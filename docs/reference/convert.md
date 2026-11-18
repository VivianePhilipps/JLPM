# Conversion

This function converts a jointLPM object to an object of type `multlcmm`
or `Jointlcmm` from lcmm package. The conversion to `multlcmm` unable
the use of the dedicated postfit and predictions functions implemented
in the lcmm package (plot, predictL, predictY, predictlink,
predictYcond, fitY, VarExpl). The conversion to `Jointlcmm` permits the
use of functions cuminc and plot (with options "baselinerisk" or
"survival").

## Usage

``` r
convert(object, to)
```

## Arguments

- object:

  an object of class `jointLPM`

- to:

  character. Either "multlcmm" or "Jointlcmm", indicating to which type
  the object should be converted.

## Value

an object of class `multlcmm` or `Jointlcmm`.
