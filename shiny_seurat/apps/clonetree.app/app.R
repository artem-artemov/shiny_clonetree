require(shiny)
require(visNetwork)
require(ape)
require(igraph)
library(Seurat)
library(ggplot2)


SR = readRDS('data/seurat.rds')

tree = read.tree('data/tree.treefile')
tree$node.label = paste0('node',1:length(tree$node.label))

leafs_select = colnames(SR)[SR@meta.data$fate=='medulla']
length(leafs_select)
tips = tree$tip.label
leafs_select_ape = paste(sapply(strsplit(leafs_select, '_'), `[`, 3), '1', sep='-')
leafs_remove = setdiff(tips, leafs_select_ape)

tree_pruned = drop.tip(tree, leafs_remove, trim.internal = FALSE)

tree_igraph <- as.igraph(tree_pruned, directed=T)

data <- toVisNetworkData(tree_igraph)

if_leafnodes <- sapply(V(tree_igraph), function(x) length(neighbors(tree_igraph,x))<=1 ) #==0
leafnodes = names(if_leafnodes)[if_leafnodes]

server <- function(input, output) {
    
  get_subtree_leafs = function(G, node, prefix = 'R30_w8.5_'){
    leafnodes_if_downstream = sapply(leafnodes, function(leaf){
        length(get.all.shortest.paths(G, V(G)[node], V(G)[leaf])$res)>0
    })
    subtree_leafs = leafnodes[leafnodes_if_downstream]
    paste0(prefix, gsub('-1', '', subtree_leafs))

  }

    
  output$network <- renderVisNetwork({
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
                           

  output$embeddingPlot = renderPlot({
      if(is.null(input$selected_node_id) | length(input$selected_node_id$nodes)==0){
          plot(
              plot(DimPlot(SR, reduction = "umap", label=T, pt.size=2))+ggtitle('clusters')+NoLegend()
          )
      }else{
          lineage = get_subtree_leafs(tree_igraph, input$selected_node_id$nodes[[1]])
          SR@meta.data$sel_lineage = ifelse(colnames(SR) %in% lineage, 1, 0)
          plot(
              FeaturePlot(SR, 'sel_lineage', pt.size=2, cols = c(rgb(.3,.3,.3,.2), rgb(1,0,0)))+ggtitle(sum(colnames(SR) %in% lineage))+NoLegend()
          )
      }
  }, height=800,width=800)
                                    #height=reactive(ifelse(!is.null(input$innerWidth),input$innerWidth*3/5,0)
        
  output$shiny_return <- renderPrint({
    input$current_node_id
  })
  output$selected_nodes <- renderPrint({
    input$selected_node_id$nodes
  })
  output$subtree_nodes <- renderPrint({
    ifelse(is.null(input$selected_node_id) | length(input$selected_node_id$nodes)==0, 'No nodes selected',  get_subtree_leafs(tree_igraph, input$selected_node_id$nodes[[1]]))
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
