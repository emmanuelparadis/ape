test_that("plot.phylo() plots node styles", {
  balanced8 <- read.tree(text = "(((a, b), (c, d)), ((e, f), (g, h)));")
  balanced10 <- read.tree(text = "((((a, b), c), (d, e)), (((f, g), h), (i, j)));")
  collapsed10 <- read.tree(text = "((((a, b), c), (d, e)), ((f, g, h), (i, j)));")
  skip_if_not_installed('vdiffr', 1.0)
  library('vdiffr')
  expect_doppelganger("node.col", function ()
    plot(balanced8, node.col = c(1, 1, 1, 1, 2, 2, 2, 3,
                                 2, 2, 1, 1, 2, 3, 3)))

  expect_doppelganger("no options", function () plot(balanced10))
  expect_doppelganger("all red", function () plot(balanced10, node.color = 'red'))
  expect_doppelganger("individual node colours", function ()
    plot(balanced10, node.color = c(rep('red', 10), rep('grey', 4),
                                    rep('red', 4), 'green')))
  expect_doppelganger("upwards", function ()
    plot(collapsed10, direction = 'upwards',
         node.color = c(rep('red', 10), rep('grey', 4), rep('red', 2), 'green', 'red'))
  )
  expect_doppelganger("nodes only", function ()
    plot(collapsed10, edge.color = 'black',
         node.color = c(rep('red', 10), rep('grey', 4),
                        rep('red', 2), 'green', 'red'))
  )

  expect_doppelganger("all options", function ()
    plot(collapsed10, direction = 'upwards',
         edge.color = hcl.colors(18, 'plasma'),
         node.color = c(rep('blue', 10), rev(hcl.colors(9, 'inferno'))),
         node.lty = c(rep('solid', 15), rep('dotted', 1), rep('dashed', 2)),
         edge.lty = 1,
         node.width = c(rep(1, 10), 1:8),
         edge.width = (1:17) / 2
         )
  )
})
