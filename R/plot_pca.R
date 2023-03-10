#' is used to plot the pairwise plots of the first 3 PCA of the gene expression (assay)
#' of a SummarizedExperiment class object.
#'
#'
#' @param pca PCA components of a SummarizedExperiment variable that will be used in the plot.
#' @param assay_names Optional selection of names of the assays to compute the PCA.
#' @param variable The variable that will be used to display and color the PCA plot
#' @param variable.name The label of the variable that will be used on the PCA plot
#' @param color The color of the variable that will be used on the PCA plot
#' @param strokeSize geom_point aesthetics
#' @param pointSize geom_point aesthetics
#' @param strokeColor geom_point aesthetics
#' @param alpha geom_point aesthetics
#'
#' @return p PCA plot of the data colored by one variable
#' @importFrom ggpubr get_legend
#' @importFrom cowplot axis_canvas ggdraw insert_xaxis_grob insert_yaxis_grob
#' @import ggplot2 scales
#' @export

plot_pca=function(
        pca,
        assay_names=NULL,
        variable,
        variable.name,
        color,
        strokeSize = .2,
        pointSize = 0.8,
        strokeColor = 'gray30',
        alpha = .5
){
    if (!is.null(assay_names)){
        normalizations=assay_names
    }else{
         normalizations=names(pca)
    }
    ppca <- lapply(
        normalizations,
        function(x){
            plot_pca_single_assay<-function(pca_x,
                                            variable,
                                            variable.name,
                                            color,
                                            strokeSize = strokeSize,
                                            pointSize = pointSize,
                                            strokeColor = strokeColor,
                                            alpha = alpha){
                pcs = pca_x$sing.val$u[,1:3]
                pc.var = pca_x$var
                pair.pcs <- utils::combn(ncol(pcs), 2)
                pList <- list()
                for(i in 1:ncol(pair.pcs)){
                    if(i == 1){
                        x <- pair.pcs[1,i]
                        y <- pair.pcs[2,i]
                        p <- ggplot(mapping = aes(
                            x = pcs[,x],
                            y = pcs[,y],
                            fill = variable)) +
                            xlab(paste0('PC', x, ' (', pc.var[x], '%)')) +
                            ylab(paste0('PC', y, ' (', pc.var[y], '%)')) +
                            geom_point(
                                aes(fill = variable),
                                pch = 21,
                                color = strokeColor,
                                stroke = strokeSize,
                                size = pointSize,
                                alpha = alpha) +
                            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            theme(
                                legend.position = "right",
                                panel.background = element_blank(),
                                axis.line = element_line(colour = "black", size = 1.1),
                                legend.background = element_blank(),
                                legend.text = element_text(size = 12),
                                legend.title = element_text(size = 14),
                                legend.key = element_blank(),
                                axis.text.x = element_text(size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.x = element_text(size = 14),
                                axis.title.y = element_text(size = 14),
                                aspect.ratio=1) +
                            guides(fill = guide_legend(override.aes = list(size = 4))) +
                            scale_fill_manual(name = variable.name, values = color)

                        le <- get_legend(p)
                    }else{
                        x <- pair.pcs[1,i]
                        y <- pair.pcs[2,i]
                        p <- ggplot(mapping = aes(
                            x = pcs[,x],
                            y = pcs[,y],
                            fill = variable)) +
                            xlab(paste0('PC', x, ' (',pc.var[x],  '%)')) +
                            ylab(paste0('PC', y, ' (',pc.var[y], '%)')) +
                            geom_point(
                                aes(fill = variable),
                                pch = 21,
                                color = strokeColor,
                                stroke = strokeSize,
                                size = pointSize,
                                alpha = alpha) +
                            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            theme(
                                panel.background = element_blank(),
                                axis.line = element_line(colour = "black", size = 1.1),
                                legend.position = "none",
                                axis.text.x = element_text(size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.x = element_text(size = 14),
                                axis.title.y = element_text(size = 14),
                                aspect.ratio=1) +
                            scale_fill_manual(values = color, name = variable.name)
                    }
                    p <- p + theme(legend.position = "none")
                    xdens <- axis_canvas(p, axis = "x")+
                        geom_density(
                            mapping = aes(
                                x = pcs[,x],
                                fill = variable),
                            alpha = 0.7,
                            size = 0.2
                        ) +
                        theme(legend.position = "none") +
                        scale_fill_manual(values = color)

                    ydens <- axis_canvas(
                        p,
                        axis = "y",
                        coord_flip = TRUE) +
                        geom_density(
                            mapping = aes(
                                x = pcs[,y],
                                fill = variable),
                            alpha = 0.7,
                            size = 0.2) +
                        theme(legend.position = "none") +
                        scale_fill_manual(name = variable.name, values = color) +
                        coord_flip()

                    p1 <- insert_xaxis_grob(
                        p,
                        xdens,
                        grid::unit(.2, "null"),
                        position = "top"
                    )
                    p2 <- insert_yaxis_grob(
                        p1,
                        ydens,
                        grid::unit(.2, "null"),
                        position = "right"
                    )
                    pList[[i]] <- ggdraw(p2)
                }
                pList[[i+1]] <- le
                return(pList)
            }
            pca_x <- pca[[x]]
            p1=plot_pca_single_assay(pca_x,variable,variable.name,color,strokeSize,pointSize,strokeColor,alpha)
            p1
        })
        names(ppca) <- normalizations

        ## Prepare plot
        p=ppca[[1]]
        if (length(normalizations)>1){
            for (n in 2:length(normalizations)){
                p=c(p,ppca[[n]])
            }
        }
        plot=do.call(
            gridExtra::grid.arrange,
            c(p,
              ncol = 4))

        return(plot=plot)
}

