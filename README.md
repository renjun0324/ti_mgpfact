

# ti_mgpfact
A Docker container for MGPfact has been developed based on [MGPfact.jl](https://github.com/renjun0324/MGPfact.jl) and [MGPfactR](https://github.com/renjun0324/MGPfactR). It is designed for one-click trajectory reconstruction and includes multiple forms of trajectory visualization.

## ti_mgpfact repository and help document

```shell
https://github.com/renjun0324/MGPfact.jl

https://github.com/renjun0324/MGPfactR

```

## pull container
```shell
docker pull renjun0324/ti_mgpfact
```

## quick start
```r

library(dynwrap)
library(dynmethods)
library(dyntoy)
library(tidyverse)
library(purrr)
library(dyno)

data("fibroblast_reprogramming_treutlein")

dataset <- wrap_expression(
  counts = fibroblast_reprogramming_treutlein$counts,
  expression = fibroblast_reprogramming_treutlein$expression
)
                               
ti_mgpfact = create_ti_method_container("renjun0324/ti_mgpfact:v0.2.0")
model = infer_trajectories(dataset_wrap, 
			    ti_mgpfact(), 
			    parameters = list(tree_method="ppt"), # or epg
			    verbose = TRUE, return_verbose = TRUE, debug = FALSE)

model$model = map(model$model, add_cell_waypoints)

modelmodel=map(modelmodel = map(modelmodel, add_cell_waypoints)
>>>>>>> d9c447cbe19d135c99638ee89897caae63513cac
metric <- map_dfr(model$model,
                  dyneval::calculate_metrics,
                  dataset = dataset,
                  metrics = c("featureimp_wcor", "him",  "F1_branches", "correlation")) 
```                        
