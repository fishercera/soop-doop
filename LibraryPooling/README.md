Library Pooling Strategy - Ralicki Method
================
Cera Fisher
5/16/2022

## Calculating pooling amounts for Illumina libraries

The purpose of this function is to allow direct pooling of a large
number of samples without first having to make 10 nM libraries for each
sample. The more libraries you have, the more time this saves. If you
have a large number of libraries *and* a wide range of final
concentrations of those libraries – we sure often do – this method will
generally allow you to keep libraries that would have been too low in
quantity to dilute to 10 nM otherwise.

For small runs, it is easier to dilute to 10nM first. Skip to the header
“Conventional Pooling Strategy” for an example of that calculation.

The original code for this was developed by Hannah Ralicki ca 2017 at
University of Connecticut, in the Jockusch Lab. Download the .Rmd file 
and the two example data files to run this for yourself in R.

### Description of input data

1.  The `quants` table

This data frame requires two columns named `Pool_ID` and
`concentration_ng_ul`. Concentration should be measured by qubit or
other fluorometry, not by Nanodrop or taken from fragment analyzer
(Bioanalyzer/TapeStation) runs. It should look like so:

| Pool_ID | concentration_ng_uL |
|:--------|:--------------------|
| A02     | 28.16666667         |
| B02     | 24.3                |
| C02     | 25.9                |
| D02     | 24.06666667         |

2.  The `regions` table

This data frame provides the information about average insert size and
also gives an estimate of how much of your library will be retained if
you do an additional size selection after pooling. A post-pooling size
selection is often desirable because it will exclude any leftover
adapter (which causes index-swapping in the Illumina ExAmp chemistry).
Additionally, because smaller fragment sizes are preferentially
sequenced by Illumina platforms, making sure that all the libraries have
the same lower-bound will help give you even sequencing coverage.

Agilent softwares for the Bioanalyzer and the TapeStation both provide
means for users to manually define regions of interest, and multiple
regions can be defined for every library. The software then
automatically calculates the average fragment size within that region,
the total molarity of the region, and what percentage of your total
library falls within that region. The user can export .csv files for all
defined regions that are roughly parallel to the data frame we use here.
For other electropherogram traces (e.g. the Agilent Fragment Analyzer),
you may have to do some work to get your information into this format.
Required fields have an asterisk.

The `regions` table has the following fields:

-   `WellID` - useful if you have multiple plates, this is the 96-well
    plate position
-   `Pool_ID` \* - The unique identifier of the library, must be the
    same as the `quants` table
-   `Plate` - the plate identifier, if you have multiple plates.
-   `SampleDescription` - A human-readable, informative name of the
    library
-   `From_bp` \*
-   `To_bp` \* - These two fields delimit the minimum and maximum length
    of a electropherogram region. You may have multiple rows for each
    library, where each row describes a different region, but you *must*
    have one row with `From_bp` = 120 and `To_bp` = 900. We use that
    range in the function to define the “total” library before size
    selection, because our experience indicates that this is the
    functional range of the sequencers we’ve used, but you can adjust
    this in the function if you like.
-   `Average_Size_bp` \* - Estimated by the fragment analysis software
    for each region
-   `Conc_pg_ul` - Estimated by the fragment analysis software, but not
    used by our function.
-   `RegionMolarity_pmol_L` - Estimated by the fragment analysis
    software, but not used by our function.
-   `percentOfTotal` \* - If you exported regions from the fragment
    analysis software, it should automatically calculate it. If not, and
    you are only using the 120-900 region range, just enter 100.
-   `region_range` \* - This is a manually created column – concatenate
    `From_bp`, -, and `To_bp` for each row.

A data-frame output from a TapeStation run with many manually defined
regions might look like this for a single sample:

| WellId | Pool_ID | Plate | SampleDescription | From_bp | To_bp | Average_size_bp | Conc_pg_uL | RegionMolarity_pmol_L | percentOfTotal | region_range |
|:-------|:--------|:------|:------------------|:--------|:------|:----------------|:-----------|:----------------------|:---------------|:-------------|
| A1     | P01-A01 | 1     | P1-A01            | 120     | 160   | 150             | 6.48       | 73.6                  | 0.55           | 120-190      |
| A1     | P01-A01 | 1     | P1-A01            | 120     | 900   | 384             | 1070       | 4650                  | 90.21          | 120-900      |
| A1     | P01-A01 | 1     | P1-A01            | 160     | 800   | 383             | 1050       | 4560                  | 89.24          | 160-800      |
| A1     | P01-A01 | 1     | P1-A01            | 250     | 800   | 394             | 993        | 4130                  | 84.04          | 250-800      |
| A1     | P01-A01 | 1     | P1-A01            | 300     | 800   | 420             | 814        | 3130                  | 68.87          | 300-800      |
| A1     | P01-A01 | 1     | P1-A01            | 350     | 800   | 454             | 596        | 2100                  | 50.47          | 350-800      |
| A1     | P01-A01 | 1     | P1-A01            | 400     | 800   | 495             | 390        | 1250                  | 33.01          | 400-800      |

There are so many rows per sample because the user defined a great many
different regions, because she was considering a variety of post-pooling
size-selection strategies. The minimum required `regions` table could
look like this for four different libraries:

**Minimum required Regions table**

| Pool_ID | From_bp | To_bp | Average_Size_bp | percentOfTotal | region_range |     |
|:--------|:--------|:------|:----------------|:---------------|:-------------|:----|
| A02     | 120     | 900   | 651             | 100            | 120-900      |     |
| B02     | 120     | 900   | 572             | 100            | 120-900      |     |
| C02     | 120     | 900   | 636             | 100            | 120-900      |     |
| D02     | 120     | 900   | 599             | 100            | 120-900      |     |

## The Library Pooling Function

First the function in its entirety:

``` r
calculatePoolingStrategy <- function(region_quants, num_libs, seq_quant, From = 120, To = 900) {
  strat <- region_quants %>% group_by(Pool_ID) %>%
    mutate(., percentLibKept=100*(percentOfTotal/percentOfTotal[region_range=="120-900"])) %>%
    filter(., From_bp == From, To_bp == To) %>%
    mutate(., strat_nM=((concentration_ng_uL*(percentLibKept/100))/(660*Average_Size_bp))*1e6) %>%
    mutate(., strat_ng_pool=((seq_quant/strat_nM)*(concentration_ng_uL*(percentLibKept/100)))/num_libs) %>%
    mutate(., strat_uL_pool = strat_ng_pool*(1/(concentration_ng_uL*(percentLibKept/100)))) %>%
    mutate(., pre_SS_ng_pooled= concentration_ng_uL*strat_uL_pool) %>%
    mutate(., n_molecules=(strat_ng_pool*6.022e23)/(Average_Size_bp*660*1e9))
  return(strat)
}
```

Then with commentary for explanation:

``` r
## This function returns a mutated region_quants table based on the number of
## libraries and the total amount of library needed for sequencing.
calculatePoolingStrategy <- function(region_quants, 
                                     num_libs, 
                                     seq_quant, 
                                     From = 120, 
                                     To = 900) {
  strat <- region_quants %>% group_by(Pool_ID) %>%
    # Add a column that calculates how much of the library is kept at each
    # pippin-prep range (uses the 120-900 range as "100%" of the library)
    mutate(., percentLibKept=100*(percentOfTotal/
                                  percentOfTotal[region_range=="120-900"])) %>%
    # Filter the region_quants table to just the entries for the region-range we
    # decided on
    filter(., From_bp == From, To_bp == To) %>%           
    # Calculate the concentration of the library we'd be keeping The
    # concentration should reflect the adapter-free library. Numerator =
    # concentration times the portion of the library kept Denominator = number
    # of picomoles based on the average fragment size
    mutate(., strat_nM=((concentration_ng_uL*(percentLibKept/100))/
                          (660*Average_Size_bp))*1e6) %>%
    # take this concentration and determine the weight needed in ng if
    # pooling to 'seq_quant' i.e. the desired nM concentration times the
    # desired volume (=the number of picomoles needed) Divide by the
    # number of samples in the pool to determine the weight per sample.
    # divide total ng needed by the number of samples to pool (n=num_libs)
    mutate(., strat_ng_pool=((seq_quant/strat_nM)*(concentration_ng_uL*
                                            (percentLibKept/100)))/num_libs) %>%
    # get the volum of each library needed for pooling
    mutate(., strat_uL_pool = strat_ng_pool*(1/(concentration_ng_uL*
                                                     (percentLibKept/100)))) %>%
    # calculate the total nanograms of each library in the pool before any size
    # selection happens
    mutate(., pre_SS_ng_pooled = concentration_ng_uL*strat_uL_pool) %>%
    # As a check, convert the post-size-selection nanograms of library in the
    # pool to the number of molecules
    mutate(., n_molecules=(strat_ng_pool*6.022e23)/(Average_Size_bp*660*1e9))
    # This number should be the same for each library -- we're shooting for
    # equimolarity.
  return(strat)
}
```

### Using the Function

#### 1. Read in quants data:

``` r
quant <- read.table(file="quants_example.txt", header=T, sep='\t')
head(quant)
```

    ##   Pool_ID concentration_ng_uL
    ## 1     A02            28.16667
    ## 2     B02            24.30000
    ## 3     C02            25.90000
    ## 4     D02            24.06667

#### 2. Read in regions data:

``` r
region_lengths <- read.table(file="regions_example.txt", header=T, sep='\t')
head(region_lengths)
```

    ##   Pool_ID From_bp To_bp Average_Size_bp percentOfTotal region_range
    ## 1     A02     120   900             651            100      120-900
    ## 2     B02     120   900             572            100      120-900
    ## 3     C02     120   900             636            100      120-900
    ## 4     D02     120   900             599            100      120-900

#### 3. Combine the data frames

``` r
region_quants <- left_join(region_lengths, quant, by = c("Pool_ID"), copy=F)    
head(region_quants)
```

    ##   Pool_ID From_bp To_bp Average_Size_bp percentOfTotal region_range
    ## 1     A02     120   900             651            100      120-900
    ## 2     B02     120   900             572            100      120-900
    ## 3     C02     120   900             636            100      120-900
    ## 4     D02     120   900             599            100      120-900
    ##   concentration_ng_uL
    ## 1            28.16667
    ## 2            24.30000
    ## 3            25.90000
    ## 4            24.06667

#### 4. Define seq_quant

We need to know the total amount of pooled library that we’re aiming
for. My sequencing center might tell me that the minimum amount of
pooled library they can take is 25 uL at 3 nM. Multiplying volume times
concentration yields amount.

``` r
seq_quant <- 25 * 3
seq_quant
```

    ## [1] 75

*I think that uL x nM resolves to femtomoles,
i.e. *m**o**l**e**s*<sup>−15</sup>, but I’m not entirely sure.* The
seq_quant gives us an amount of total library that we’re going to divide
across the number of libraries we have. You can make this amount much
larger if you have plenty of library and don’t want to try to pipette
small amounts. `seq_quant` is the variable we play with the most when
working out pooling strategies.

#### 5. Run the library pooling function

``` r
pooling_strat <- calculatePoolingStrategy(region_quants = region_quants, 
                                          num_libs = 4, # Well, that's all I have!  
                                          seq_quant = seq_quant, 
                                          From = 120, # 120 and 900 are defaults, 
                                          To = 900)   # but here I set them explicitly

pooling_strat
```

    ## # A tibble: 4 x 13
    ## # Groups:   Pool_ID [4]
    ##   Pool_ID From_bp To_bp Average_Size_bp percentOfTotal region_range
    ##   <chr>     <int> <int>           <int>          <int> <chr>       
    ## 1 A02         120   900             651            100 120-900     
    ## 2 B02         120   900             572            100 120-900     
    ## 3 C02         120   900             636            100 120-900     
    ## 4 D02         120   900             599            100 120-900     
    ## # ... with 7 more variables: concentration_ng_uL <dbl>, percentLibKept <dbl>,
    ## #   strat_nM <dbl>, strat_ng_pool <dbl>, strat_uL_pool <dbl>,
    ## #   pre_SS_ng_pooled <dbl>, n_molecules <dbl>

And let’s look at what this strategy is suggesting that I do with my
physical hands in the lab:

``` r
### How many uL of each library do I add to my pool? 

pooling_strat[,c(1,11)]
```

    ## # A tibble: 4 x 2
    ## # Groups:   Pool_ID [4]
    ##   Pool_ID strat_uL_pool
    ##   <chr>           <dbl>
    ## 1 A02             0.286
    ## 2 B02             0.291
    ## 3 C02             0.304
    ## 4 D02             0.308

I’m obviously not going to be able to pipette such small amounts! This
is why we adjust the `seq_quant` value the most when we’re figuring out
how to pool. These particular libraries are all extremely uniform, but
when libraries have wider ranges of concentration and/or sizing, we
often find it most fruitful to separate them out into different groups,
pool them separately, and then pool the pools.

I’ll just re-run this again with a more reasonable `seq_quant`:

``` r
pooling_strat <- calculatePoolingStrategy(region_quants = region_quants, 
                                          num_libs = 4, # Well, that's all I have!  
                                          seq_quant = 500, 
                                          From = 120, # 120 and 900 are defaults, 
                                          To = 900)   # but here I set them explicitly

pooling_strat[,c(1,11)]
```

    ## # A tibble: 4 x 2
    ## # Groups:   Pool_ID [4]
    ##   Pool_ID strat_uL_pool
    ##   <chr>           <dbl>
    ## 1 A02              1.91
    ## 2 B02              1.94
    ## 3 C02              2.03
    ## 4 D02              2.05

Those are perfectly comfortable amounts to pipette. But wait! These
volumes will not add up to my minimum 25 uL, not even close! That’s
because the libraries are much more concentrated than 3 nM. In fact, we
have a good estimate of their current, native nM:

``` r
pooling_strat$strat_nM
```

    ## [1] 65.55571 64.36745 61.70192 60.87587

My `seq_quant` of 500 could be thought of as 50 uL of 10 nM library,
which is totally comfortable. You can pipette the calculated amounts,
and then make up the volume to 50 uL with 1 M Tris. It won’t be spot on
10 nM, but it will be about right.

Again, the value of this direct pooling strategy really shines when you
have to pool across multiple plates of libraries, some of which may not
have given a final yield \> 10 nM. With just a few well-behaved
libraries, one probably is better off figuring out how to dilute to 10
nM and then just add equal amounts of the 10 nM library. As a kindness
to the reader, who is expected to primarily be myself, I will do that
here.

### Conventional Pooling Strategy

I have already calculated the native nM of my libraries as described
above, but let’s recalculate it to prove to ourselves we understand it.

I’ll re-use the region_quants dataframe:

``` r
head(region_quants)
```

    ##   Pool_ID From_bp To_bp Average_Size_bp percentOfTotal region_range
    ## 1     A02     120   900             651            100      120-900
    ## 2     B02     120   900             572            100      120-900
    ## 3     C02     120   900             636            100      120-900
    ## 4     D02     120   900             599            100      120-900
    ##   concentration_ng_uL
    ## 1            28.16667
    ## 2            24.30000
    ## 3            25.90000
    ## 4            24.06667

Then use the molecular weight of DNA to get nanomolarity:

``` r
region_quants <- region_quants %>% 
    mutate(., pre_dilution_nM=((concentration_ng_uL)/
                          (660*Average_Size_bp))*1e6) 
```

Let’s say that I want to make 10 uL of 10 nM library for each of these
libraries:

``` r
vol_wanted <- 10
nM_wanted <- 10 

# V1 x C1 = V2 x C2
# V1 = (V2 x C2) / C1

region_quants <- region_quants %>% 
  mutate(., vol_wanted = vol_wanted, nM_wanted = nM_wanted) %>%
  mutate(., vol_needed = (vol_wanted * nM_wanted)/pre_dilution_nM) %>%
  mutate(., vol_buffer = vol_wanted - vol_needed)

data.frame(Pool_ID = region_quants$Pool_ID, vol_uL = region_quants$vol_needed, vol_buffer_uL = region_quants$vol_buffer)
```

    ##   Pool_ID   vol_uL vol_buffer_uL
    ## 1     A02 1.525420      8.474580
    ## 2     B02 1.553580      8.446420
    ## 3     C02 1.620695      8.379305
    ## 4     D02 1.642687      8.357313

This is also a completely comfortable pipetting amount, and results in
having more than enough 10 nM library to pool for a minimum of 25 uL.
