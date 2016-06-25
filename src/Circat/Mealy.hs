{-# LANGUAGE CPP #-}
{-# LANGUAGE ExistentialQuantification, TupleSections, TypeOperators #-}
{-# LANGUAGE ScopedTypeVariables, Arrows #-}
{-# LANGUAGE ConstraintKinds, TypeFamilies #-}
{-# LANGUAGE Rank2Types #-} -- for tests
{-# OPTIONS_GHC -Wall #-}

-- {-# OPTIONS_GHC -fno-warn-unused-imports #-} -- TEMP
-- {-# OPTIONS_GHC -fno-warn-unused-binds   #-} -- TEMP

#define CircuitConstraint

-- #define Semantics

----------------------------------------------------------------------
-- |
-- Module      :  Circat.Mealy
-- Copyright   :  (c) 2014 Tabula, Inc.
-- 
-- Maintainer  :  conal@tabula.com
-- Stability   :  experimental
-- 
-- Mealy machines
----------------------------------------------------------------------

module Circat.Mealy
  ( Mealy(..),runMealy,ArrowCircuit(..),ArrowCircuitT
  , scanl, scan, iterateU, foldSP
  , scanl', scan', iterateU'
  , sumS, productS
  , sumSP, productSP, dotSP
  , sumPS, productPS, dotPS
  , mac, sumMac
#ifdef Semantics
  , asStreamFun, asArrow
#endif
--   , mealyC  -- experiment
  , asFun, mealy
  , MStream, scanlStr, scanStr
  , double, maybeM
  ) where

import Prelude hiding (id,(.),sum,product,scanl)
import Control.Category
import Control.Arrow
import Control.Applicative (liftA2)
import Data.Monoid ((<>),Sum(..),Product(..))
import Data.Foldable (Foldable(..),sum,product)
import Data.Tuple (swap)
import Data.List (unfoldr)
import GHC.Exts (Constraint)

import Circat.Misc (dup) -- ,Unit,(:*)
import qualified Circat.Misc as M

#ifdef CircuitConstraint
import Circat.Circuit (GenBuses)
#endif

#ifdef Semantics
import Data.Stream (Stream(..))
import qualified Control.Arrow.Operations as Op
import Control.Arrow.Transformer.Stream
#endif

{--------------------------------------------------------------------
    Experiment with constraining s
--------------------------------------------------------------------}

type GS a = (GenBuses a, Show a)

#ifdef CircuitConstraint

-- type C = GenBuses
type C a = GS a

-- The Show constraint is helpful for debugging. Perhaps remove later.

#else

class C s

instance C ()
instance (C a, C b) => C (a,b)
instance C Bool
instance C Int

#endif

--------------------------------------------------------------------

data Mealy a b = forall s. C s => Mealy ((a,s) -> (b,s)) s

-- We could probably generalize Mealy to an arrow transformer.

instance Category Mealy where
  id = Mealy id ()
  Mealy g t0 . Mealy f s0 = Mealy h (s0,t0)
   where
     h (a,(s,t)) = (c,(s',t'))
      where
        (b,s') = f (a,s)
        (c,t') = g (b,t)
  {-# INLINE id #-}
  {-# INLINE (.) #-}

instance Arrow Mealy where
  arr f = Mealy (first f) ()
  Mealy f s0 *** Mealy g t0 = Mealy h (s0,t0)
   where
     h ((a,b),(s,t)) = ((c,d),(s',t'))
      where
        (c,s') = f (a,s)
        (d,t') = g (b,t)
  first  f = f *** id
  second g = id *** g
  {-# INLINE arr #-}
  {-# INLINE (***) #-}
  {-# INLINE first #-}
  {-# INLINE second #-}

instance ArrowChoice Mealy where
  Mealy f s0 +++ Mealy g t0 = Mealy h (s0,t0)
   where
     h (Left  a,(s,t)) = (Left  c,(s',t)) where (c,s') = f (a,s)
     h (Right b,(s,t)) = (Right d,(s,t')) where (d,t') = g (b,t)
  left  f = f +++ id
  right g = id +++ g
  {-# INLINE (+++) #-}
  {-# INLINE left #-}
  {-# INLINE right #-}

instance ArrowLoop Mealy where
  loop (Mealy f s0) = Mealy (loop (rot . f . rot)) s0
   where
     rot :: ((x,y),z) -> ((x,z),y)
     rot ((x,y),z) = ((x,z),y)

class ArrowLoop k => ArrowCircuit k where
  type CircuitKon k a :: Constraint
  type CircuitKon k a = ()
  delay :: CircuitKon k a => a -> k a a

instance ArrowCircuit Mealy where
  type CircuitKon Mealy a = C a
  delay = Mealy swap
  {-# INLINE delay #-}

-- delay a = Mealy (\ (new,old) -> (old,new)) a

type ArrowCircuitT k a = (ArrowCircuit k, CircuitKon k a)

-- A sort of semantic function.
runMealy :: Mealy a b -> [a] -> [b]
runMealy (Mealy f s0) = go s0
 where
   go _ []     = []
   go s (a:as) = b : go s' as where (b,s') = f (a,s)

#ifdef Semantics

type StreamFun = StreamArrow (->)

asStreamFun :: Mealy a b -> StreamFun a b
asStreamFun (Mealy f s0) = StreamArrow (go s0)
 where
   go s (Cons a as) = Cons b (go s' as) where (b,s') = f (a,s)

-- Much more generally:

asArrow :: Op.ArrowCircuit k =>
           Mealy a b -> (a `k` b)
asArrow (Mealy f s0) = loop (arr f . second (Op.delay s0))

#endif

-- | Impossible interpretation in '(->)' (which doesn't really have delay) but
-- possible in other categories. To be inlined and re-interpreted.
asFun :: Mealy a b -> (a -> b)
asFun (Mealy h s0) = mealy h s0
-- asFun (Mealy h s0) = M.loop (h . second (M.delay s0))
{-# INLINE asFun #-}

mealy :: ((a,s) -> (b,s)) -> s -> (a -> b)
mealy h s0 = M.loop (h . second (M.delay s0))
{-# INLINE mealy #-}

{--------------------------------------------------------------------
    Some other standard instances for arrows
--------------------------------------------------------------------}

instance Functor (Mealy a) where
  fmap = flip (>>^)
  {-# INLINE fmap #-}

instance Applicative (Mealy a) where
  pure b = arr (const b)
  fs <*> xs = uncurry ($) <$> (fs &&& xs)
  {-# INLINE pure #-}
  {-# INLINE (<*>) #-}

instance Monoid b => Monoid (Mealy a b) where
  mempty = pure mempty
  mappend = liftA2 mappend

-- | Scan via Mealy

scanl :: C b => (b -> a -> b) -> b -> Mealy a b
scanl op = Mealy (\ (a,b) -> dup (b `op` a))
{-# INLINE scanl #-}

-- scanl op = Mealy (dup . uncurry (flip op))

-- TODO: Generalize scanl to ArrowCircuit

scan :: (Monoid m, C m) => Mealy m m
scan = scanl (<>) mempty
{-# INLINE scan #-}

iterateU :: C a => (a -> a) -> a -> Mealy () a
iterateU f = scanl (\ a () -> f a)
{-# INLINE iterateU #-}

-- | Like 'scanl' but yields the initial value first.
-- Maybe replace 'scanl' by this variant.
scanl' :: C b => (b -> a -> b) -> b -> Mealy a b
scanl' op = Mealy (\ (a,b) -> (b, b `op` a))
{-# INLINE scanl' #-}

scan' :: (Monoid m, C m) => Mealy m m
scan' = scanl' (<>) mempty
{-# INLINE scan' #-}

iterateU' :: C a => (a -> a) -> a -> Mealy () a
iterateU' f = scanl' (\ a () -> f a)
{-# INLINE iterateU' #-}

sumS :: (Num a, GS a) => Mealy a a
sumS = arr getSum . scan . arr Sum
{-# INLINE sumS #-}

-- sumS = Mealy (\ (a,tot) -> dup (tot+a)) 0
-- sumS = Mealy (dup . uncurry (+)) 0

productS :: (Num a, GS a) => Mealy a a
productS = arr getProduct . scan . arr Product
{-# INLINE productS #-}

-- productS = Mealy (\ (a,tot) -> dup (tot*a)) 0
-- productS = Mealy (dup . uncurry (*)) 0

-- Sequential-of-parallel

foldSP :: (Foldable f, Monoid m, C m) => Mealy (f m) m
foldSP = scan . arr fold
{-# INLINE foldSP #-}

-- Equivalently,

-- foldSP = Mealy (dup . uncurry mappend) mempty . arr fold
-- foldSP = Mealy (dup . uncurry mappend) mempty . Mealy (first fold) ()
-- foldSP = Mealy (\ (as,tot) -> dup (fold as <> tot)) mempty
-- foldSP = scanl (\ as tot -> fold as <> tot) mempty

sumSP :: (Functor f, Foldable f, Num a, C a) => Mealy (f a) a
sumSP = arr getSum . foldSP . arr (fmap Sum)
{-# INLINE sumSP #-}

productSP :: (Functor f, Foldable f, Num a, C a) => Mealy (f a) a
productSP = arr getProduct . foldSP . arr (fmap Product)
{-# INLINE productSP #-}

-- | Dot product with sum over o (outer) and product over i (inner). Note that
-- this sum/product choice is a bit odd, but I think it's a good fit for
-- sequential-of-parallel inner products. Otherwise, one can simply transpose
-- first with @'dotSP' . 'arr' 'transpose'@.
dotSP :: (Foldable o, Functor o, Foldable i, Num a, C a) =>
         Mealy (o (i a)) a
dotSP = sumSP . arr (fmap product)
{-# INLINE dotSP #-}

_dotSP2 :: (Foldable o, Functor o, Foldable i, Num a, C a) =>
           Mealy (o (i a)) a
_dotSP2 = Mealy (\ (pas,tot) -> dup (dot pas + tot)) 0
 where
   dot = sum . fmap product

-- Parallel-of-sequential

-- foldPS :: (Foldable f, Monoid (f a), Monoid a, GS (f a)) => Mealy (f a) a
-- foldPS = arr fold . scan
-- {-# INLINE foldPS #-}

sumPS :: (Foldable f, Num (f a), Num a, GS (f a)) =>
         Mealy (f a) a
sumPS = arr sum . sumS
-- sumPS = Mealy (\ (t,tot) -> let tot' = tot+t in (sum tot',tot')) 0
{-# INLINE sumPS #-}

productPS :: (Foldable f, Num (f a), Num a, GS (f a)) =>
             Mealy (f a) a
productPS = arr product . productS
-- productPS = Mealy (\ (t,tot) -> let tot' = tot+t in (product tot',tot')) 0
{-# INLINE productPS #-}

dotPS :: (Foldable o, Functor o, Foldable i, Num (o a), Num a, GS (o a)) =>
         Mealy (o (i a)) a
dotPS = sumPS . arr (fmap product)
{-# INLINE dotPS #-}

mac :: (Foldable f, Num a, GS a) =>
       Mealy (f a) a
mac = sumS . arr product
-- mac = Mealy (\ (as,tot) -> dup (tot+product as)) 0
{-# INLINE mac #-}

sumMac :: (Foldable o, Foldable i, Num (i a), Num a, GS (i a)) =>
          Mealy (o (i a)) a
sumMac = arr sum . mac
-- sumMac = Mealy (\ (as,tot) -> let tot' = tot+product as in (sum tot',tot')) 0
{-# INLINE sumMac #-}

{--------------------------------------------------------------------
    An encoding of synchronous streams
--------------------------------------------------------------------}

type MStream = Mealy ()

-- type StreamFun a b = MStream a -> MStream b

scanlStr :: C b => (b -> a -> b) -> b -> MStream a -> MStream b
scanlStr op e ma = scanl op e . ma

scanStr :: (Monoid m, C m) => MStream m -> MStream m
scanStr = scanlStr (<>) mempty

{--------------------------------------------------------------------
    Misc experiments
--------------------------------------------------------------------}

double :: Mealy a b -> Mealy (a,a) (b,b)
double (Mealy g s0) = Mealy h s0
 where
   h ((a,a'),s) = ((b,b'),s'')
    where
      (b ,s' ) = g (a ,s )
      (b',s'') = g (a',s')

-- double :: Mealy a b -> Mealy a b
-- double (Mealy f s0) = Mealy h (s0,t0,False)
--  where
--    h (a,(s,t,i)) = ...
--     where
--       ... = f s

maybeM :: Mealy a b -> Mealy (Maybe a) (Maybe b)
maybeM (Mealy g s0) = Mealy h s0
 where
   h (Nothing,s) = (Nothing,s)
   h (Just a ,s) = first Just (g (a,s))

{--------------------------------------------------------------------
    Examples
--------------------------------------------------------------------}

-- Counter

counter0 :: (C a, Num a) => Mealy () a
counter0 = Mealy (\ ((),n) -> (n,n+1)) 0

counter1 :: (ArrowCircuitT k a, Num a) => () `k` a
counter1 = loop (arr (\ ((),n) -> (n,n+1)) . second (delay 0))

-- Also works when the composed functions are swapped.

counter2 :: (ArrowCircuitT k a, Num a) => () `k` a
counter2 = proc () -> do rec old <- delay 0 -< new
                             new <- returnA -< old + 1
                         returnA -< old

testCounter :: Mealy () Int -> [Int]
testCounter counter = runMealy counter (replicate 10 ())

_testCounters :: [[Int]]
_testCounters = testCounter <$> [counter0,counter1,counter2]

-- *Circat.Mealy> _testCounters
-- [[0,1,2,3,4,5,6,7,8,9],[0,1,2,3,4,5,6,7,8,9],[0,1,2,3,4,5,6,7,8,9]]

adder0 :: (C a, Num a) => Mealy a a
adder0 = Mealy (\ (old,a) -> dup (old+a)) 0

adder1 :: (ArrowCircuitT k a, Num a) => a `k` a
adder1 = loop (arr (\ (a,tot) -> dup (tot+a)) . second (delay 0))

-- With arrow notation. Oops. Exclusive.
adder2 :: (ArrowCircuitT k a, Num a) => a `k` a
adder2 = proc a -> do rec tot <- delay 0 -< tot + a
                      returnA -< tot

-- With arrow notation. Inclusive.
adder3 :: (ArrowCircuitT k a, Num a) => a `k` a
adder3 = proc a -> do rec old <- delay 0 -< new
                          new <- returnA -< old + a
                      returnA -< new

-- Using let. Inclusive.
adder4 :: (ArrowCircuitT k a, Num a) => a `k` a
adder4 = proc a -> do rec old <- delay 0 -< new
                          let new = old + a
                      returnA -< new

testAdder :: Mealy Int Int -> [Int]
testAdder adder = runMealy adder [1..10]

_testAdders :: [[Int]]
_testAdders = testAdder <$> [adder0,adder1,adder2,adder3,adder4]

fibL1 :: Num a => [a]
fibL1 = fib' (0,1)
 where
   fib' (a,b) = a : fib' (b,a+b)

fibL1' :: Num a => [a]
fibL1' = unfoldr (\ (a,b) -> Just (a,(b,a+b))) (0,1)

fibL2 :: Num a => [a]
fibL2 = 0 : 1 : zipWith (+) fibL2 (tail fibL2)

fib0 :: (C a, Num a) => Mealy () a
fib0 = Mealy (\ ((),(a,b)) -> (a,(b,a+b))) (0,1)

fib1 :: (C a, Num a) => Mealy () a
fib1 = loop (arr (\ ((),(a,b)) -> (a,(b,a+b))) . second (delay (0,1)))

fib2 :: (C a, Num a) => Mealy () a
fib2 = proc () -> do rec (a,b) <- delay (0,1) -< (b,a+b)
                     returnA -< a

fib3 :: (C a, Num a) => Mealy () a
fib3 = proc () -> do rec a <- delay 0 -< b
                         b <- delay 1 -< a+b
                     returnA -< a

testFib :: Mealy () a -> [a]
testFib fib = runMealy fib (replicate 10 ())


_testFibs :: [[Int]]
_testFibs = (take 10 <$> [fibL1,fibL1',fibL2]) ++ (testFib <$> [fib0,fib1,fib2,fib3])
