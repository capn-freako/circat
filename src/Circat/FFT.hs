{-# LANGUAGE GADTs #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}

----------------------------------------------------------------------
-- |
-- Module      :  Circat.FFT
-- Copyright   :  (c) 2015 David Banas
-- License     :  BSD3
-- 
-- Maintainer  :  capn.freako@gmail.com
-- Stability   :  experimental
-- 
-- Fast Fourier Transform (FFT), as a class, via Circat
----------------------------------------------------------------------

module Circat.FFT where

import Prelude hiding ({- id,(.), -}foldl,foldr,sum,product,zipWith,reverse,and,or,scanl,minimum,maximum)

import Control.Applicative
import Control.Arrow
import Data.Complex
import Data.Foldable (Foldable, sum, foldl')
import TypeUnary.Nat (IsNat, Nat(..), nat, N2, N3, N4, N5)  -- , N6)

import Circat.Scan (lproducts, LScan)
import qualified Circat.Pair as P
import qualified Circat.RTree as RT
import Circat.RTree (bottomSplit)

type RTree = RT.Tree

-- | FFT, as a class
-- (The LScan constraint comes from the use of 'lproducts', in 'addPhase'.)
class (LScan f) => FFT f a where
    fft  :: f a -> f a  -- ^ Computes the FFT of a functor.

-- Note that this definition of the FFT instance for Pair assumes DIT.
-- How can we eliminate this assumption and make this more general?
instance (RealFloat a, Applicative f, Foldable f, Num (f (Complex a)), FFT f (Complex a)) => FFT P.Pair (f (Complex a)) where
    fft = P.inP (uncurry (+) &&& uncurry (-)) . P.secondP addPhase . fmap fft

instance (IsNat n, RealFloat a) => FFT (RTree n) (Complex a) where
    fft = fft' nat
        where   fft' :: (RealFloat a) => Nat n -> RTree n (Complex a) -> RTree n (Complex a)
                fft' Zero     = id
                fft' (Succ _) = inDIT fft
                    where   inDIT g  = RT.toB . g . bottomSplit

-- | Adds the proper phase adjustments to a functor containing Complex RealFloats,
-- and instancing Num.
addPhase :: (Applicative f, Foldable f, LScan f, RealFloat a, Num (f (Complex a))) => f (Complex a) -> f (Complex a)
addPhase = liftA2 (*) id phasor
  where phasor f = fst $ lproducts (pure phaseDelta)
          where phaseDelta = cis ((-pi) / fromIntegral n)
                n          = flen f

-- | Gives the "length" (i.e. - number of elements in) of a Foldable.
-- (Soon, to be provided by the Foldable class, as "length".)
flen :: (Foldable f) => f a -> Int
flen = foldl' (flip ((+) . const 1)) 0

