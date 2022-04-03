require(shiny)
require(visNetwork)
require(ape)
require(igraph)
library(Seurat)
library(ggplot2)

#SR = readRDS('/mnt/data/rd')

SR = readRDS('~/Adrenal/adrenal/data/Seurat/adrenal.human.seurat.rds')
SR = SR[,SR@meta.data$orig.ident=='R30_w85']
gc()
SR = SR[,SR@meta.data$orig.ident=='R30_w85']
gc()
#cellnames = paste(sapply(strsplit(colnames(SR), '_'), `[`, 3), '1', sep='-')
#head(cellnames)


tree = read.tree('/home/artem/Adrenal/mito_tree/mito_genotypes.pos2000cells100.fasta.treefile')
tree$node.label = paste0('node',1:length(tree$node.label))

leafs_select = colnames(SR)[SR@meta.data$fate=='medulla']
length(leafs_select)
tips = tree$tip.label
leafs_select_ape = paste(sapply(strsplit(leafs_select, '_'), `[`, 3), '1', sep='-')
leafs_remove = setdiff(tips, leafs_select_ape)

tree_pruned = drop.tip(tree, leafs_remove, trim.internal = FALSE)
#plot(tree_pruned)

#tree_igraph <- as.igraph(tree, directed=T)
tree_igraph <- as.igraph(tree_pruned, directed=T)

data <- toVisNetworkData(tree_igraph)

if_leafnodes <- sapply(V(tree_igraph), function(x) length(neighbors(tree_igraph,x))<=1 ) #==0
leafnodes = names(if_leafnodes)[if_leafnodes]

server <- function(input, output) {
    
  get_subtree_leafs = function(G, node, prefix = 'R30_w8.5_'){
    #if_leafnodes <- sapply(V(G), function(x) length(neighbors(G,x))<=1 ) #==0
    #leafnodes = names(if_leafnodes)[if_leafnodes]
    leafnodes_if_downstream = sapply(leafnodes, function(leaf){
        length(get.all.shortest.paths(G, V(G)[node], V(G)[leaf])$res)>0
    })
    subtree_leafs = leafnodes[leafnodes_if_downstream]
    paste0(prefix, gsub('-1', '', subtree_leafs))

  }

    
  output$network <- renderVisNetwork({
    # minimal example
    #nodes <- data.frame(id = 1:3)
    #edges <- data.frame(from = c(1,2), to = c(1,3))
    
    #visNetwork(nodes, edges)
      
    #tree_str = "((A,(B,(X,(Y,Z)))),C);"
      #"((A:1,(B:1,(X,Y):1):2):3,C:5);"
      #'((A:1,B:2):3,C:5);'

    #tree <- read.tree(text = tree_str)   #library(ape)
    #tree = read.tree('/home/artem/Adrenal/mito_tree/mito_genotypes.pos2000cells100.fasta.treefile')
    #tree$node.label = paste0('node',1:length(tree$node.label))
    #tree_igraph <- as.igraph(tree)         #library(igraph)

    #data <- toVisNetworkData(tree_igraph)
    visNetwork(nodes = data$nodes, edges = data$edges, height = "500px", maxNodeSize = 0.1) %>% 
          visEdges(arrows = "from") %>% 
          visNodes(size=10) %>%
          visHierarchicalLayout() %>%
              visPhysics(stabilization = FALSE) %>%
              visEdges(smooth = FALSE) %>%
          #visIgraphLayout(layout = "layout.reingold.tilford", circular=T) %>% #layout_as_tree
          visInteraction(hover = TRUE) %>%
          visEvents(hoverNode = "function(nodes) {
              Shiny.onInputChange('current_node_id', nodes);
          ;}") %>%
          visEvents(select = "function(nodes) {
              Shiny.onInputChange('selected_node_id', nodes);
          ;}")

  })
                           
  #output$distPlot <- renderPlot({
  #  hist(rnorm(100), col = 'darkgray', border = 'white')
  #})

  #output$subtree_leafs = get_subtree_leafs(tree_igraph, input$selected_node_id$nodes[[1]])
  output$embeddingPlot = renderPlot({
      if(is.null(input$selected_node_id)){
          plot(
              plot(DimPlot(SR, reduction = "umap", label=T))+ggtitle('clusters')+NoLegend()
          )
      }else{
          lineage = get_subtree_leafs(tree_igraph, input$selected_node_id$nodes[[1]])
          #output$subtree_leafs
          #print(lineage)
          SR@meta.data$sel_lineage = ifelse(colnames(SR) %in% lineage, 1, 0)
          plot(
              FeaturePlot(SR, 'sel_lineage', pt.size=0.5, cols = c(rgb(.3,.3,.3,.2), rgb(1,0,0)))+ggtitle(sum(colnames(SR) %in% lineage))+NoLegend()
          )
      }
  }, height=800,width=800)
                                    #height=reactive(ifelse(!is.null(input$innerWidth),input$innerWidth*3/5,0)
        
  #output$distPlot <- renderPlot({
  #  hist(rnorm(100), col = 'darkgray', border = 'white')
  #})
  output$shiny_return <- renderPrint({
    input$current_node_id
  })
  output$selected_nodes <- renderPrint({
    input$selected_node_id$nodes
  })
  output$subtree_nodes <- renderPrint({
    get_subtree_leafs(tree_igraph, input$selected_node_id$nodes[[1]])
  })

}

ui <- fluidPage(
    titlePanel("Clonal tree"),

    fluidRow(

    column(4,
        visNetworkOutput("network"),
        wellPanel(
            verbatimTextOutput("shiny_return"),
            verbatimTextOutput("selected_nodes"),
            verbatimTextOutput("subtree_nodes")
            #sliderInput("obs", "Number of observations:",  
            #          min = 1, max = 1000, value = 500)
        )       
    ),

    column(8,
      plotOutput("embeddingPlot")
    )
    )
)

shinyApp(ui = ui, server = server)
