library(gt)

# Convert data to data frame
mat <- as.data.frame(t(lres))

df_num_cols <- sapply(mat, is.numeric)
mat[, df_num_cols] <- lapply(mat[, df_num_cols], function(x) sprintf("%.2f", x))
mat[mat == "NA"] <- ""
mat[mat == "0.00"] <- ""
mat[mat == "-0.00"] <- ""

colnames(mat) <- c("CESD", "SWL")

# Add row names as a new column
mat <- cbind(RowNames = rownames(mat), mat)

# Add a new column for row groups
mat$RowGroup <- c(rep("Affect shifts", 10), rep("Mean/sd", 4))

# Create the table
tbl <- gt(mat, rowname_col = "RowNames", groupname_col = "RowGroup")

tbl <- cols_label(tbl,
    CESD = "CESD", SWL = "SWL"
)
# tbl <- fmt_number(tbl, decimals=2)
tbl <- tab_header(tbl, title = md("**LASSO coefficients**"))
# tbl <- tab_caption(tbl, md("**Table 1**. LASSO coefficients"))
tbl <- tab_options(tbl,
    row.striping.include_table_body = TRUE,
    row_group.as_column = TRUE,
    table.border.top.style = "solid",
    table.border.top.width = px(1),
    table.border.bottom.style = "solid",
    table.border.bottom.width = px(1)
)

# Group columns under "Study 1" and "Study 2"
# tbl <- tab_spanner(tbl, label = "Study 1", columns = 2:7)
# tbl <- tab_spanner(tbl, label = "Study 2", columns = 8:12)
