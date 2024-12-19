# load data for tests
spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
islet_poly <- reconstructShapeDensityImage(spe,
    marks = "cell_category",
    image_col = "image_name", image_id = "E04", mark_select = "islet", dim = 500
)
metrics_matrix <- totalShapeMetrics(islet_poly)


test_that("shapeMetrics computes correct shape metrics", {
    # Create a simple polygon
    matrix_R <- matrix(c(
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 0, 0, 0,
        0, 1, 1, 0, 0, 1, 1, 0, 0,
        0, 1, 1, 0, 0, 1, 1, 0, 0,
        0, 1, 1, 1, 1, 1, 0, 0, 0,
        0, 1, 1, 0, 1, 1, 0, 0, 0,
        0, 1, 1, 0, 0, 1, 1, 0, 0,
        0, 1, 1, 0, 0, 1, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 9, byrow = TRUE)

    poly_R <- binaryImageToSF(matrix_R, xmin = 0, xmax = 1, ymin = 0, ymax = 1)

    # Test function
    metrics <- shapeMetrics(poly_R)

    expect_type(metrics, "list")
    expect_named(metrics, c(
        "Area", "Compactness", "Eccentricity", "Circularity",
        "Solidity", "Curl", "fibreLength", "fibreWidth"
    ))
    expect_gt(metrics$Area, 0)
    expect_gt(metrics$Compactness, 0)

    expect_gte(metrics$Eccentricity, 0)
    expect_lte(metrics$Eccentricity, 1)

    expect_gte(metrics$Circularity, 0)
    expect_lte(metrics$Circularity, 1)

    expect_gte(metrics$Solidity, 0)
    expect_lte(metrics$Solidity, 1)

    expect_gt(metrics$Curl, 0)
})

test_that("totalShapeMetrics computes shape metrics matrix correctly", {
    # Test function
    expect_type(metrics_matrix, "double")
    expect_true(ncol(metrics_matrix) > 0)
    expect_true(nrow(metrics_matrix) >= 1)
})

test_that("meanShapeMetrics computes mean shape metrics correctly", {
    # Test function
    mean_metrics <- meanShapeMetrics(metrics_matrix)
    expect_true(ncol(mean_metrics) == 1)
})

### Invalid inputs

test_that("shapeMetrics throws an error for invalid input", {
    invalid_input <- list("not a polygon")

    expect_error(shapeMetrics(invalid_input), "'sfPoly' must be a valid sfc object")
})

test_that("totalShapeMetrics throws an error for invalid input", {
    invalid_input <- list("not a multipolygon")

    expect_error(totalShapeMetrics(invalid_input), "'sfInput' must be a valid sf object")
})

test_that("meanShapeMetrics throws an error for invalid input", {
    invalid_input <- list("not a matrix")

    expect_error(meanShapeMetrics(invalid_input), "'totalShapeMetricMatrix' must be a matrix")
})

