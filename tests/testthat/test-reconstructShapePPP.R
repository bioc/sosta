# Create a dummy point pattern
ppp <- spatstat.geom::ppp(x = runif(100),
                          y = runif(100),
                          window = spatstat.geom::owin(c(0, 1), c(0, 1)))

# load SPE object
spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)

# Test the reconstruction function
islet_poly <- reconstructShapeDensityImage(
    spe,
    marks = "cell_category",
    image_col = "image_name",
    image_id = "E04",
    mark_select = "islet",
    dim = 500
)

# Test the reconstruction on SPE object
spe_sel <- spe[, spe[["image_name"]] %in% c("E02", "E03", "E04")]
all_islets <- reconstructShapeDensitySPE(
    spe_sel,
    marks = "cell_category",
    image_col = "image_name",
    mark_select = "islet",
    bndw = sigma,
    thres = 0.0025
)


test_that("reconstructShapeDensity returns valid polygons", {
    # Reconstruct polygons with valid parameters
    result <- reconstructShapeDensity(ppp, mark_select = NULL, dim = 100)

    expect_s3_class(result, "sf")
    expect_true(any(st_geometry_type(result) == "POLYGON"))
    expect_false(any(st_is_empty(result)))
})

test_that("reconstructShapeDensity handles invalid thresholds", {
    # Test low threshold
    expect_error(reconstructShapeDensity(ppp, thres = 0, dim = 500),
                 "Threshold too low")

    # Test high threshold
    expect_error(reconstructShapeDensity(ppp, thres = 1E5, dim = 500),
                 "Threshold too high")
})

test_that("shapeIntensityImage generates valid plots", {
    # Test the intensity image function
    p <- shapeIntensityImage(
        spe,
        marks = "cell_category",
        image_col = "image_name",
        image_id = "E04",
        mark_select = "islet"
    )

    expect_s3_class(p, "ggplot")
})

test_that("reconstructShapeDensityImage returns polygons from SpatialExperiment",
          {
              expect_s3_class(islet_poly, "sf")
              expect_true(any(st_geometry_type(islet_poly) == "POLYGON"))
              expect_false(any(st_is_empty(islet_poly)))
          })

test_that("reconstructShapeDensitySPE handles multiple images", {
    expect_s3_class(all_islets, "sf")
    expect_true(any(st_geometry_type(all_islets) == "POLYGON"))
    expect_false(any(st_is_empty(all_islets)))
})

test_that("estimateReconstructionParametersSPE estimates valid parameters",
          {
              # Test the estimation function
              res <- estimateReconstructionParametersSPE(
                  spe = spe,
                  marks = "cell_category",
                  image_col = "image_name",
                  mark_select = "islet",
                  nimages = 2,
                  dim = 500,
                  plot_hist = FALSE
              )

              expect_true(is.data.frame(res))
              expect_named(res, c("img", "bndw", "thres"))
              expect_true(all(res$bndw > 0))
              expect_true(all(res$thres > 0))
          })

test_that("reconstructShapeDensity handles invalid input types", {
    # Invalid input (non-ppp object)
    expect_error(
        reconstructShapeDensity("not_a_ppp"),
        "'ppp' must be an object of class 'ppp'"
    )

    expect_error(reconstructShapeDensity(ppp, dim = -500),
                 "'dim' must be a single, positive, numeric value")
})

test_that("estimateReconstructionParametersSPE handles edge cases", {

    # Test with fewer images than requested
    expect_error(
        estimateReconstructionParametersSPE(
            spe = spe_sel,
            marks = "cell_category",
            image_col = "image_name",
            mark_select = "islet",
            nimages = 10,
            dim = 500,
            plot_hist = FALSE
        ),
        "must be smaller or equal to the number of images"
    )
})
