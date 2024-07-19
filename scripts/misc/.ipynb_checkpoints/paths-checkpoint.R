# Variables
ver <- "07-05-2024"

## PATHS and common variables
version_dir <- glue::glue("{ver}")

proj_dir <- paste0(here::here(),"/")

## Create directory for plots
plt_dir <- glue::glue("{ver}/plots_{ver}")

## Create directory for RDS objects
robj_dir <- glue::glue("{ver}/R_objects_{ver}")
markers_dir <- glue::glue("{ver}/Markers_{ver}")
