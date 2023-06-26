# plotting_utils
## Contains custom plotting themes and wrapper functions...


## **scThemes**
This wrapper function takes in a few parameters for your figures (font sizes, etc.) and returns a list of ggplot themes for a few plot types commonly used in single-cell analysis/figures.

#### Example usage:
```
source("/path/to/DWM_utils/plotting_utils/scThemes.R")

scTheme <- scThemes(
  small.font = 6,
  big.font = 8,
  line.width = 1,
  pt.size = 0.1
)

ggplot(
  df,
  aes(
    x=x,
    y=y,
    color=color
  )
)+
geom_point() +
scThemes$scatter
```

## **McKolors**
Just a compiled list of color palettes that I use. Includes palettes from R and python packages, plus some other random ones. Sources below.
#### Color palette sources & resources
- `scico` - [[link](https://github.com/thomasp85/scico)]
- `seaborn` - [[link](https://seaborn.pydata.org/tutorial/color_palettes.html)]
- `pals` - [[link](https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html)]
- `ggsci` -  [[link](https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html)]
- `MetBrewer` - [[link](https://github.com/BlakeRMills/MetBrewer)]
- `wesanderson` - [[link](https://github.com/karthik/wesanderson)]
- "gadenbuie/ggpomological"
- "m-clark/NineteenEightyR"
- "oompaBase"

https://coolors.co/palettes/popular

#### To-Add:
- National Parks Palettes: https://github.com/kevinsblake/NatParksPalettes
- MetBrewer - [link](https://github.com/BlakeRMills/MetBrewer)
- "gadenbuie/ggpomological"

#### Example usage (in R):
```
mckolors <- read.csv("/path/to/DWM_utils/plotting_utils/McKolors_v1.csv") %>%
  as.list() %>%
  lapply(
    FUN=function(X) X[X!=""]
  )

  ggplot(
    df,
    aes(
      x=x,
      y=y,
      color=color
    )
  ) +
  geom_point() +
  scThemes$scatter +
  scale_color_manual(
    values=mckolors$ldw29
  )
```
