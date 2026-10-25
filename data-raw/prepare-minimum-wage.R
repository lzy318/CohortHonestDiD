# Reproduces data/minimum_wage.rda. The input is the county panel of the
# minimum-wage application in Callaway and Sant'Anna (2021), as used in
# Liu (2026), "Cohort-Anchored HonestDiD". The input file itself is not
# included in this repository: save it as an R data file at
# data-raw/min_wage_CS.rds (or change source_path below) and run this script
# from the package root. The script adds convenience names; it does not change
# the source data.
source_path <- file.path("data-raw", "min_wage_CS.rds")
if (!file.exists(source_path)) stop("Study input not found: ", source_path)
minimum_wage <- readRDS(source_path)
minimum_wage$id <- minimum_wage$countyreal
minimum_wage$Y <- minimum_wage$lemp
minimum_wage$D <- as.integer(minimum_wage$first.treat > 0 &
                             minimum_wage$year >= minimum_wage$first.treat)
minimum_wage$state <- as.integer(minimum_wage$FIPS %/% 1000)
minimum_wage$cohort <- minimum_wage$first.treat
rownames(minimum_wage) <- NULL
dir.create("data", showWarnings = FALSE)
save(minimum_wage, file = file.path("data", "minimum_wage.rda"), compress = "xz")
