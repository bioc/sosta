test_that("st_feature_axes computes major and minor axis lengths correctly", {
    # Create a simple polygon
    matrix_R <- matrix(c(
        0, 1, 1, 0,
        0, 1, 1, 0,
        0, 1, 1, 0,
        0, 1, 1, 0
    ), nrow = 4, byrow = TRUE)

    poly_R <- binaryImageToSF(matrix_R, xmin = 0, xmax = 1, ymin = 0, ymax = 1)

    # Test function
    axes <- st_feature_axes(poly_R)

    expect_type(axes, "list")
    expect_named(axes, c("majorAxisLength", "minorAxisLength"))
    expect_equal(axes$majorAxisLength, 1)
    expect_equal(axes$minorAxisLength, 0.5)
})

test_that("st_calculateCurvature computes curvature metrics correctly", {
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
    curvature <- st_calculateCurvature(poly_R)

    expect_type(curvature, "list")
    expect_named(curvature, c("meanAbsCurv", "totalAbsCurv", "minCurv", "maxCurv"))
    expect_true(curvature$meanAbsCurv > 0)
    expect_true(curvature$totalAbsCurv > 0)
})

test_that("st_calculateShapeCurl computes curl correctly", {
    # Create a simple polygon
    matrix_R <- matrix(c(
        1, 1, 1, 1, 1, 0,
        1, 1, 0, 0, 1, 1,
        1, 1, 0, 0, 1, 1,
        1, 1, 1, 1, 1, 0,
        1, 1, 0, 1, 1, 0,
        1, 1, 0, 0, 1, 1,
        1, 1, 0, 0, 1, 1
    ), nrow = 7, byrow = TRUE)

    poly_R <- binaryImageToSF(matrix_R, xmin = 0, xmax = 1, ymin = 0, ymax = 1)

    # Test function
    curl <- st_calculateShapeCurl(poly_R)

    expect_type(curl, "list")
    expect_named(curl, c("Curl", "fibreLength", "fibreWidth"))
    expect_true(curl$Curl >= 0)
    expect_true(curl$fibreLength > 0)
    expect_true(curl$fibreWidth > 0)
})

### Invalid inputs

test_that("st_feature_axes throws an error for invalid input", {
    invalid_input <- list("not a polygon")

    expect_error(st_feature_axes(invalid_input), "'sfPoly' must be a valid sfc object")
})

test_that("st_calculateCurvature throws an error for invalid input", {
    invalid_input <- list("not a polygon")

    expect_error(st_calculateCurvature(invalid_input), "'sfPoly' must be a valid sfc object")
})

test_that("st_calculateShapeCurl throws an error for invalid input", {
    invalid_input <- list("not a polygon")

    expect_error(st_calculateShapeCurl(invalid_input), "'sfPoly' must be a valid sfc object")
})
