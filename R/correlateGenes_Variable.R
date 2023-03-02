#' is used to find
#'
#'
#' @param se.obj PCs of the dataset that will be used in the plot
#' @param assay All the available assays for the data (i.e. normalizations methods)
#' @param variable The regression variable that will be computed to the PCAof the data
#' @param apply.log The label of the variable that will be used on the PCA plot
#' @param a The number of components of the PCA used to compute the regression
#' @param rho The number of components of the PCA used to compute the regression
#' @param check.top.genes The number of components of the PCA used to compute the regression
#' @param top.genes.no The number of components of the PCA used to compute the regression
#'
#' @return list List containing the association plot and the computed regression

#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assay
#' @importFrom Rfast correls transpose
#' @export

correlateGenes_Variable <- function(
    se.obj,
    assay,
    variable,
    apply.log = TRUE,
    method = "pearson",
    a = 0.05,
    rho = 0,
    check.top.genes = TRUE,
    top.genes.no = 3
)
{
    ### check the inputs
    if( anyNA(se.obj@colData[, variable]) )
        message('There is NA in the variable, please remove them and re-run the function')
    if( anyNA(se.obj@colData[, variable]) )
        message('There is NA in the assay data, please remove them and re-run the function')
    ### log transformation
    if(apply.log){
        y <- log2(assay(x = se.obj, assay) + 1)
        message('Performing log transformation on the data + 1')
    }
    else{
        y <- assay(x = se.obj, assay)
    }
    ### correlation
    message(paste0(
        'Performing ' ,
        method,
        ' correlation between individual genes and the ',
        variable,
        ' variable')
    )
    corr.genes.var <- correls(
        y = se.obj@colData[ , variable],
        x = transpose(y),
        type = method,
        a = a ,
        rho = rho
    )
    row.names(corr.genes.var) <- row.names(se.obj)
    if(check.top.genes){
        corr.genes.var <- corr.genes.var[
            order(
                corr.genes.var[ , 'correlation'],
                decreasing = TRUE,
                na.last = TRUE) , ]
        ### positive correlation
        p.high <- as.data.frame(t(y[row.names(corr.genes.var)[c(1:top.genes.no)], ]))
        p.high$var = se.obj@colData[ , variable]
        p.high <-  pivot_longer(-var, names_to = 'genes', values_to = 'expr')
        p.high <- ggplot(p.high, aes(x = var, y = expr)) +
            geom_point() +
            ylab('Gene expression (log)') +
            xlab(variable) +
            facet_wrap(~genes) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', size = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                strip.text.x = element_text(size = 10),
                plot.title = element_text(size = 16)
            )
        ### negative correlation
        p.low <- as.data.frame(t(y[row.names(corr.genes.var)[c(c(nrow(corr.genes.var)-c(top.genes.no -1)): nrow(corr.genes.var))], ]))
        p.low$var <- se.obj@colData[ , variable]
        p.low <- pivot_longer(-var, names_to = 'genes', values_to = 'expr')
        p.low <- ggplot(p.low, aes(x = var, y = expr)) +
            geom_point() +
            ylab('Gene expression (log)') +
            xlab(variable) +
            facet_wrap(~genes) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', size = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                strip.text.x = element_text(size = 10),
                plot.title = element_text(size = 16)
            )
        return(list(
            corr.genes.var = corr.genes.var[row.names(se.obj), ],
            plot.low = p.low,
            plot.high = p.high))

    }
    else{
        return(corr.genes.var = corr.genes.var[row.names(se.obj), ])
    }
}
