#' Draw aggregation plot.
#' 
#' PlotAPA
#' @description Draw aggregation plot from aggregation matrices.
#' @param aggregatedMtx <matrix>: The aggregated matrix.
#' @param trim <numeric>: A number between 0 and 100 thaht give the percentage of triming in matrices.
#' @param colMin <matrix>: The minimal value in color scale. If Null automaticaly find.
#' @param colMid <matrix>: The middle value in color scale. If Null automaticaly find.
#' @param colMax <matrix>: The mximal value in color scale. If Null automaticaly find.
#' @param colMinCond <matrix>: Avalaible for plotting differantial aggregation. The minimal value in color scale in the classsical aggregation plot. If Null automaticaly find.
#' @param colMaxCond <matrix>: Avalaible for plotting differantial aggregation. The maxiaml value in color scale in the classsical aggregation plot. If Null automaticaly find.
#' @return None
#' @examples
#' library(GenomicED)
#' data(aggreg.mtx)
#' 
#' 
#' PlotAPA(
#'     aggregatedMtx                  = aggreg.mtx,
#'     trim             = 20,
#'     colMin          = -2,
#'     colMid               = 0,
#'     colMax          = 2,
#'     colMinCond = 0,
#'     colMaxCond = 2
#' )
#' 
#' 

    # aggregatedMtx = NULL, title = NULL,
    # trim = 0, tails = "both",
    # colMin = NULL, colMid = NULL,
    # colMax = NULL, colBreaks = NULL,
    # blurPass = 0, boxKernel = NULL,
    # kernSize = NULL,
    # stdev = 0.5, loTri = NULL,
    # colors = NULL,
    # na.value = "#F2F2F2",
    # colorScale = "linear",
    # bias = 1, paletteLength = 51

PlotAPA = function(aggregatedMtx = NULL, trim=0, colMin=NULL, colMid=NULL, colMax=NULL, colMinCond=NULL, colMaxCond=NULL){
    .ggDensity <- function(data.lst=NULL, colour.col=NULL, mean.bln=TRUE, title=NULL){
        data.lst_tbl <- lapply(seq_along(data.lst),function(element.ndx){
            return(tibble::tibble(
                value = data.lst[[element.ndx]],
                class = factor(names(data.lst)[[element.ndx]])
            ))
            }) 
        data.tbl <- dplyr::bind_rows(data.lst_tbl)
        if(is.null(colour.col)){
            colour.col <- Hue(length(data.lst)) |> stats::setNames(names(data.lst))
        }
        plot.gp <- ggplot2::ggplot(data.tbl, ggplot2::aes(x=data.tbl$value, fill=data.tbl$class, colour=data.tbl$class)) + 
            ggplot2::geom_density(alpha=0.1) +
            ggplot2::scale_color_manual(values = colour.col)+
            ggplot2::scale_fill_manual(values = colour.col) +
            ggplot2::labs(title=title) 
        if (mean.bln){
            data.tbl <- dplyr::group_by(data.tbl, class = data.tbl$class)
            mu.tbl <-  dplyr::summarise(data.tbl, grp.mean = mean(data.tbl$value))
            plot.gp <- plot.gp + 
                ggplot2::geom_vline(data = mu.tbl, ggplot2::aes(xintercept = mu.tbl$grp.mean, colour = mu.tbl$class), linetype = "dashed")
        }
        return(plot.gp)
    }
    # Differential or not?
        differential.bln <- !is.null(attributes(aggregatedMtx)$matrices)
        if(differential.bln){
            colors = NULL
        }else{
            colors = viridis(255)
        }
    # Plot
        # Auto Scale
            plot.gp <- ggAPA(
                aggregatedMtx=aggregatedMtx, 
                colors=colors,
                title=ifelse(differential.bln,
                    yes="Agregation of differential matrices",
                    no="Agregation")
            ) + ggplot2::labs(subtitle="scale (auto), center()")
            print(plot.gp)
        # Auto Scale + Center
            if(!is.null(colMid)){
                plot.gp <- ggAPA(
                    aggregatedMtx=aggregatedMtx, 
                    colors=colors,
                    colMid=colMid,
                    title=ifelse(differential.bln,
                        yes="Agregation of differential matrices",
                        no="Agregation")
                ) + ggplot2::labs(subtitle=paste0("scale (auto), center(",colMid,")"))
                print(plot.gp)
            }
        # Trim Scale + Center
            if(!is.null(trim) && 0<trim){
                plot.gp <- ggAPA(
                    aggregatedMtx=aggregatedMtx, 
                    colors=colors,
                    trim=trim,
                    colMid=colMid,
                    title=ifelse(differential.bln,
                        yes="Agregation of differential matrices",
                        no="Agregation")
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trim,"%), center(",colMid,")"))
                print(plot.gp)
            }
        # MinMax Scale + Center
            if(!is.null(colMin) || !is.null(colMax)){
                plot.gp <- ggAPA(
                    aggregatedMtx=aggregatedMtx, 
                    colors=colors,
                    colMin=colMin,
                    colMid=colMid,
                    colMax=colMax,
                    title=ifelse(differential.bln,
                        yes="Agregation of differential matrices",
                        no="Agregation")
                ) + ggplot2::labs(subtitle=paste0("scale (",colMin,";",colMax,"), center(",colMid,")"))
                print(plot.gp)
            }
    if (differential.bln){
        # Pval + Auto Scale
            if(!is.null(attributes(aggregatedMtx)$matrices$pVal) && sum(!is.na(attributes(aggregatedMtx)$matrices$pVal))>=3){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$pVal,
                    colors=YlOrRd(9),
                    title = "-log10(p.values)"
                ) + ggplot2::labs(subtitle="scale (auto), center()")
            }else{
                plot.gp <- ggplot2::ggplot() +
                    ggplot2::theme_void() +
                    ggplot2::annotate("text", x = 1, y = 1,
                        label = "Not enough pval computed to plot a pval matrix (<3) or nothing significant")
            }
            print(plot.gp)
        # FiltPval + Auto Scale + Center
            if(!is.null(attributes(aggregatedMtx)$matrices$aggDiffPvalFilt) && sum(!is.na(attributes(aggregatedMtx)$matrices$pVal))>=3){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$aggDiffPvalFilt,
                    colors=colors,
                    colMid=colMid,
                    title = paste0("Agregation of differential matrices")
                ) + ggplot2::labs(subtitle=paste0("filtred by p.values, scale (auto), center(",colMid,")"))
            }else{
                plot.gp <- ggplot2::ggplot() +
                    ggplot2::theme_void() +
                    ggplot2::annotate("text", x = 1, y = 1,
                    label = "Not enough pval computed to plot a pval matrix (<3) or nothing significant")
            }
            print(plot.gp)
        # FiltPval + Trim Scale + Center
            if(!is.null(attributes(aggregatedMtx)$matrices$aggDiffPvalFilt) && sum(!is.na(attributes(aggregatedMtx)$matrices$pVal))>=3){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$aggDiffPvalFilt,
                    colors=colors,
                    trim=trim,
                    colMid=colMid,
                    title = paste0("Agregation of differential matrices")
                ) + ggplot2::labs(subtitle=paste0("filtred by p.values, scale (rm ",trim,"%), center(",colMid,")"))
            }else{
                plot.gp <- ggplot2::ggplot() +
                    ggplot2::theme_void() +
                    ggplot2::annotate("text", x = 1, y = 1,
                    label = "Not enough pval computed to plot a pval matrix (<3) or nothing significant")
            }
            print(plot.gp)
        # Delta + Auto Scale + Center
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggDelta, 
                colors=colors,
                colMid=colMid,
                title="Differential of agregated matrices"
            ) + ggplot2::labs(subtitle=paste0("scale (auto), center(",colMid,")"))
            print(plot.gp)
        # Delta + Trim Scale + Center
            if(!is.null(trim) && 0<trim){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$aggDelta, 
                    colors=colors,
                    trim=trim,
                    colMid=colMid,
                    title="Differential of agregated matrices"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trim,"%), center(",colMid,")"))
                print(plot.gp)
            }
        # Delta + Auto Scale + Center
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrectedDelta, 
                colors=colors,
                colMid=colMid,
                title="Differential of corrected agregated matrices"
            ) + ggplot2::labs(subtitle=paste0("scale (auto), center(",colMid,")"))
            print(plot.gp)
        # Delta + Trim Scale + Center
            if(!is.null(trim) && 0<trim){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrectedDelta, 
                    colors=colors,
                    trim=trim,
                    colMid=colMid,
                    title="Differential of corrected agregated matrices"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trim,"%), center(",colMid,")"))
                print(plot.gp)
            }
        # Control + Auto Scale
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
                colors=viridis(51),
                title="Agregation control"
            ) + ggplot2::labs(subtitle="scale (auto)")
            print(plot.gp)
        # Control + Trim Scale
            if(!is.null(trim) && 0<trim){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
                    colors=viridis(51),
                    trim=trim,
                    title="Agregation control"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trim,"%)"))
                print(plot.gp)
            }
        # Control + MinMax Scale
            if(!is.null(colMinCond) || !is.null(colMaxCond)){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
                    colors=viridis(51),
                    colMin=colMinCond,
                    colMax=colMaxCond,
                    title="Agregation control"
                ) + ggplot2::labs(subtitle=paste0("scale (",colMinCond,";",colMaxCond,")"))
                print(plot.gp)
            }
        # Condition + Auto Scale
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$agg, 
                colors=viridis(51),
                title="Agregation"
            ) + ggplot2::labs(subtitle="scale (auto)")
            print(plot.gp)
        # Condition + Trim Scale
            if(!is.null(trim) && 0<trim){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$agg, 
                    colors=viridis(51),
                    trim=trim,
                    title="Agregation"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trim,"%)"))
                print(plot.gp)
            }
        # Condition + MinMax Scale
            if(!is.null(colMinCond) || !is.null(colMaxCond)){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$agg, 
                    colors=viridis(51),
                    colMin=colMinCond,
                    colMax=colMaxCond,
                    title="Agregation"
                ) + ggplot2::labs(subtitle=paste0("scale (",colMinCond,";",colMaxCond,")"))
                print(plot.gp)
            }
        # Corrected condition + Auto Scale
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrected, 
                colors=viridis(51),
                title="Agregation corrected"
            ) + ggplot2::labs(subtitle="scale (auto)")
            print(plot.gp)
        # Corrected condition + Trim Scale
            if(!is.null(trim) && 0<trim){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrected, 
                    colors=viridis(51),
                    trim=trim,
                    title="Agregation corrected"
                ) + ggplot2::labs(subtitle=paste0("scale (rm ",trim,"%)"))
                print(plot.gp)
            }
        # Corrected condition + MinMax Scale
            if(!is.null(colMinCond) || !is.null(colMaxCond)){
                plot.gp <- ggAPA(
                    aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrected, 
                    colors=viridis(51),
                    colMin=colMinCond,
                    colMax=colMaxCond,
                    title="Agregation corrected"
                ) + ggplot2::labs(subtitle=paste0("scale (",colMinCond,";",colMaxCond,")"))
                print(plot.gp)
            }
        # Grouped Scale(Condition & Control)
            colBreaks.num <- BreakVector(
                x.num=c(c(attributes(aggregatedMtx)$matrices$agg), c(attributes(aggregatedMtx)$matrices$aggCtrl)),
                n.num=51)
        # Control + Grouped Scale(with Condition)
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
                colors=viridis(51),
                colBreaks.num=colBreaks.num,
                title="Agregation control"
            ) + ggplot2::labs(subtitle="scale (grouped with condition)")
            print(plot.gp)
        # Condition + Grouped Scale(with Control)
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$agg, 
                colors=viridis(51),
                colBreaks.num=colBreaks.num,
                title="Agregation"
            ) + ggplot2::labs(subtitle="scale (grouped with control)")
            print(plot.gp)
        # Grouped Scale(Corrected condition & Control)
            colBreaks.num <- BreakVector(
                x.num=c(c(attributes(aggregatedMtx)$matrices$aggCorrected), c(attributes(aggregatedMtx)$matrices$aggCtrl)),
                n.num=51)
        # Control + Grouped Scale(with corrected condition)
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
                colors=viridis(51),
                colBreaks.num=colBreaks.num,
                title="Agregation control"
            ) + ggplot2::labs(subtitle="scale (grouped with condition corrected)")
            print(plot.gp)
        # Corrected condition + Grouped Scale(with Control)
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrected, 
                colors=viridis(51),
                colBreaks.num=colBreaks.num,
                title="Agregation corrected"
            ) + ggplot2::labs(subtitle="scale (grouped with control)")
            print(plot.gp)
        # Density

            plot.gp <- .ggDensity(
                data.lst=list(differential=stats::na.omit(c(aggregatedMtx))),
                colour.col=Hue(3)[[1]],
                title="Agregation of differential matrices density")
            print(plot.gp)

            data.lst <- list(
                control=stats::na.omit(c(attributes(aggregatedMtx)$matrices$aggCtrl)),
                condition=stats::na.omit(c(attributes(aggregatedMtx)$matrices$agg))
            )
            plot.gp <- .ggDensity(
                data.lst=data.lst,
                colour.col=Hue(3)[2:3],
                title="Condition and control densities")
            print(plot.gp)
            
            plot.gp <- .ggDensity(
                data.lst=list(deltaCorrected=stats::na.omit(c(attributes(aggregatedMtx)$matrices$aggCorrectedDelta))),
                colour.col=Hue(3)[[1]],
                title="Differential of corrected agregated matrices density")
            print(plot.gp)

            data.lst <- list(
                control=stats::na.omit(c(attributes(aggregatedMtx)$matrices$aggCtrl)),
                correctedCondition=stats::na.omit(c(attributes(aggregatedMtx)$matrices$aggCorrected))
            )
            plot.gp <- .ggDensity(
                data.lst=data.lst,
                colour.col=Hue(3)[2:3],
                title="Condition corrected and control densities")
            print(plot.gp)
    }
    # Attributes
            grid::grid.newpage()
            attr.ndx <- aggregatedMtx |>
                attributes() |>
                names() |>
                NotIn(c("dim","matrices","interactions", "dimnames")) # NotIn to replace by %ni%
            attr.lst <- attributes(aggregatedMtx)[attr.ndx]
            attr.lst$aggregationMethod <- function(pxl){pxl[is.na(pxl)]<-0;mean(pxl,na.rm=TRUE)}
            attr1.ndx <- attr.lst |>
                lapply(class) |>
                unlist() |>
                NotIn(c("matrix", "list","GInteractions","function")) # NotIn to replace by %ni% 
            attr1.lst <- attr.lst[attr1.ndx] |>
                        lapply(as.character) |>
                        unlist()
            attr2.ndx <- unlist(lapply(attr.lst, class)) %in% "function"
            attr2.lst <- attr.lst[attr2.ndx] |>
                        lapply(function(function.fun){
                            function.chr <- deparse(function.fun)
                            function.chr[3:(length(function.chr)-1)] |>
                            paste0(collapse=";\n")
                            return(gsub(" ","",function.chr))
                            }) |>
                        unlist()
            attr.lst <- c(attr1.lst, attr2.lst)
            attr.tbl  <- tibble::as_tibble_col(attr.lst) |>
                tibble::add_column(name=names(attr.lst)) |>
                tibble::column_to_rownames("name")
            gridExtra::grid.table(attr.tbl)
}