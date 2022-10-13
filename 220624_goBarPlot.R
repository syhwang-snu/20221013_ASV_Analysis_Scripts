


goBarPlot <- function(x, plot.title='',pvalue = 'p.adjust' ,showCategory = 12){
    
    category <- min(nrow(x), showCategory)
    df <- x %>% arrange(p.adjust) %>% 
        slice_min(p.adjust, n = category, with_ties = FALSE)  %>% 
        mutate(Description = factor(Description, level = rev(Description)), 
               log10padj = log10(pvalue),
               position = rep(1,category))
    
    p1 <- df %>% 
        ggplot(aes(x = Description ,y= Count)) +
        geom_col(fill = 'transparent', color = 'black', size = 0.7, width = 0.8) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 35)) +
        coord_flip() +
        theme_void() + 
        theme(
            axis.title.x = element_text(size = 12, face = "bold") ,
            axis.title.y = element_blank() ,
            axis.text = element_text(face = "bold", size = 12),
            axis.text.x = element_text(face = "bold", size = 11),
            axis.text.y = element_text(hjust = 0.95),
            legend.position = 'none',
            plot.background = element_rect(fill = 'transparent', color = NA),
            panel.background = element_rect(fill = 'transparent',color = 'black'),
            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm") ,
            legend.background = element_rect(fill='transparent',color = NA),
            legend.box.background = element_rect(fill='transparent', color = NA)
        ) 
    
    p2 <- df %>% 
        ggplot(aes(x = Description )) +
        geom_text(aes(y = position, label = round(log10padj,2)), hjust = 1.5,vjust = -0, size = 4.5, fontface = 'bold') + coord_flip() + theme_void()+ geom_vline(xintercept = seq(0.5,(category + 1.5), 1)) +
        labs(y = '', x= NULL) + ggtitle('log10(P-Value)') + theme(text = element_text(face = 'bold'), 
                                                                  plot.title = element_text(hjust = 0.5))
    
    p3 <- df %>% 
        ggplot(aes(x = Description )) +
        geom_text(aes(y = position, label = str_wrap(geneID_FC_top_genes, width = 100))) + 
        geom_vline(xintercept = seq(0.5,(category + 1.5), 1)) +
        coord_flip() + theme_void()+
        labs(y = '', x= NULL) + ggtitle('Genes') + theme(text = element_text(face = 'bold'), 
                                                         plot.title = element_text(hjust = 0.5))
    
    p <- p1 + p2 + p3 + plot_layout(widths = c(1,1,5)) + 
        patchwork::plot_annotation(title = plot.title, theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.4),
                                                            plot.background = element_rect(fill = 'transparent', color = NA),
                                                            panel.background = element_rect(fill = 'transparent',color = NA),
                                                            legend.background = element_rect(fill='transparent', color = NA),
                                                            legend.box.background = element_rect(fill='transparent', color = NA)
        ))
    return(p)

}


