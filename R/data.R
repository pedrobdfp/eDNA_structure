#' Example eDNA metabarcoding dataset
#'
#' @description
#' A simulated eDNA metabarcoding dataset representing a coastal fisheries
#' survey. The dataset has known ground-truth community structure (K = 3
#' communities) driven by depth and latitude, making it ideal for learning
#' the package and verifying that your installation works correctly.
#'
#' Load it via [get_example_data()].
#'
#' @format A named list with three elements:
#' \describe{
#'   \item{`counts`}{An integer matrix (60 rows × 25 columns).
#'     Rows are sampling stations (named `STN_001` through `STN_060`).
#'     Columns are taxa (25 species of NE Pacific marine vertebrates and
#'     cephalopods). Values are simulated read counts.}
#'   \item{`covariates`}{A data frame with 60 rows and 5 columns:
#'     \describe{
#'       \item{`sample_id`}{Station identifier.}
#'       \item{`latitude`}{Decimal degrees N (range: 39.1–41.4°N).}
#'       \item{`depth`}{Sample depth in meters (50, 150, or 300 m).}
#'       \item{`year`}{Sampling year (2019, 2021, or 2023).}
#'       \item{`depth_bin`}{Depth stratum: `"shallow"` (≤100 m) or `"deep"` (>100 m).}
#'     }
#'   }
#'   \item{`true_params`}{A list containing the true (data-generating) parameters:
#'     `K = 3` communities, `pi` (3 × 25 composition matrix), and `gamma`
#'     (60 × 3 membership probability matrix). These are provided so you can
#'     compare model estimates to the truth.}
#' }
#'
#' @section Data-generating process:
#' The data were generated from a Dirichlet-Multinomial mixture model with
#' K = 3 communities:
#' \itemize{
#'   \item **Community 1** (shallow + southern): dominated by small pelagic fish
#'     (*Engraulis mordax*, *Sardinops sagax*).
#'   \item **Community 2** (deep + northern): dominated by groundfish and
#'     mesopelagics (*Merluccius productus*, *Sebastes* spp., *Anoplopoma fimbria*).
#'   \item **Community 3** (transitional, reference): mixed composition including
#'     mesopelagic lanternfish (*Stenobrachius leucopsarus*) and squid.
#' }
#' Community membership is driven by depth (strong effect, β ≈ ±2) and
#' latitude (moderate effect, β ≈ ±1.5). Overdispersion alpha = 3.
#' Read depth per sample ranges from 2,186 to 7,938.
#'
#' @source Simulated data; taxon names are real NE Pacific species but counts
#'   are synthetic.
#'
#' @examples
#' data <- get_example_data()
#' dim(data$counts)          # 60 x 25
#' head(data$covariates)
#'
#' @name example_edna
#' @docType data
#' @keywords datasets
NULL
