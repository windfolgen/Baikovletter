# Baikovletter
compute symbol letters for multi-loop planar or non-planar MPL Feynman integral families

## contents

### BaikovAll.wl
a Mathematica package which generates Baikov representations for a given Feynman integral family by recursively integrating Baikov variables.
### BaikovLetter.wl
a Mathematica package which generates rational letters and algebraic letters for a given multiloop planar MPL Feynman integral family.
### UsageofBaikovLetter.nb
a Mathematica notebook used to demonstrate the usage of BaikovLetter.wl
### UsageofBaikovAll.wl
a Mathematica notebook used to demonstrate the usage of BaikovAll.wl
### CDE/
CDE matrices and UT basis for two nontrivial examples: two-loop two-mass pentabox and three-loop two-massive-leg ladder

## Description
This package __BaikovLetter.wl__ contains two major part:

First, it will generate all possible leading singularities in Baikov representation which will be the candidates for rational letters of an integral family. In this step, it will also pick out all possible leading singularities which are square roots. These will serve as constraints for the next step to determine all algebraic letters. In this step we will use the third-part program __Singular__ to identify all spurious letters. You can download it in the following site. For Mac users, it can be easily installed by [Homebrew](https://brew.sh).

> <https://www.singular.uni-kl.de/index.html>

Second, it will generate all algebraic letters for this integral family using the method described by our paper [arXiv: 2401.07632](http://arxiv.org/abs/2401.07632).

At last, we put canonical differential equations and UT bases of two nontrivial integrals families in the file folder CDE/.
