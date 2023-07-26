library(ggplot2)
library(ggpubr)
load("plot_data.rdata")

ui = navbarPage("TCGADEPMAP", tabPanel("Gene essentiality across lineages",fluidPage(
  titlePanel("Gene essentiality across lineages"),
  # Create a new Row in the UI for selectInputs
  fluidRow(
    column(4,
           selectInput("data",
                       "Data:",
                       c("DEPMAP",
                         unique(c("TCGADEPMAP(exp.only)","GTEXDEPMAP(exp.only)","PDXEDEPMAP(exp.only)","TCGADEPMAP(integrated)"))))),
    column(4,
           mainPanel(uiOutput("esssel")))
  ),
  mainPanel(
    plotOutput(outputId = "distPlot"),
    uiOutput("down")
    
  ))),
  tabPanel("Gene essentiality vs mut/del/amp",fluidPage(
    titlePanel("Gene essentiality vs mut/del/amp"),
    radioButtons("muttype", "mut/del/amp type:",
                 c("Mutation" = "mut",
                   "Deletion" = "del",
                   "Amplification" = "amp")),
    # Create a new Row in the UI for selectInputs
    fluidRow(
      column(4,
             selectInput("mutdata",
                         "Data:",
                         c("DEPMAP",
                           unique(c("TCGADEPMAP(exp.only)","PDXEDEPMAP(exp.only)","TCGADEPMAP(integrated)"))))
      ),
      column(4,
             mainPanel(uiOutput("esssel2"))),
      column(4,
             mainPanel(uiOutput("mutsel")))
    ),
    mainPanel(
      plotOutput(outputId = "mutPlot"),
      uiOutput("down1")
    )))
)

server = function(input, output) {
  xdata = reactive({
    x_data = c("")
    if(input$gene != "")
    {
      x_gene = input$gene
      if(input$data == "DEPMAP")
      {
        x_data = depmap_plot[,c(x_gene,"lineage")]
      }
      else if(input$data == "TCGADEPMAP(exp.only)")
      {
        x_data = tcga_plot[,c(x_gene,"lineage")]
      }
      else if(input$data == "GTEXDEPMAP(exp.only)")
      {
        x_data= gtex_plot[,c(x_gene,"lineage")]
      }
      else if(input$data == "PDXEDEPMAP(exp.only)")
      {
        x_data = pdx_plot[,c(x_gene,"lineage")]
      }
      else if(input$data == "TCGADEPMAP(integrated)")
      {
        x_data = tcga_multi_plot[,c(x_gene,"lineage")]
      }
    }
    x_data
  })
  # Filter data based on selections
  output$esssel = renderUI({
    sel_list = c("")
    if(input$data == "DEPMAP")
    {
      sel_list = colnames(pdx_plot)[!grepl("^MUT[.]",colnames(pdx_plot))]
    }
    else if(input$data == "TCGADEPMAP(exp.only)")
    {
      sel_list = colnames(pdx_plot)[!grepl("^MUT[.]",colnames(pdx_plot))]
    }
    else if(input$data == "PDXEDEPMAP(exp.only)")
    {
      sel_list = colnames(pdx_plot)[!grepl("^MUT[.]",colnames(pdx_plot))]
    }
    else if(input$data == "GTEXDEPMAP(exp.only)")
    {
      sel_list = colnames(pdx_plot)[!grepl("^MUT[.]",colnames(pdx_plot))]
    }
    else if(input$data == "TCGADEPMAP(integrated)")
    {
      #sel_list = colnames(tcga_multi_plot)[!(grepl("^MUT[.]",colnames(tcga_multi_plot))|grepl("^DEL[.]",colnames(tcga_multi_plot))|grepl("^AMP[.]",colnames(tcga_multi_plot)))]
      sel_list = colnames(tcga_multi_plot)
    }
    sel_list = setdiff(sel_list,"lineage")
    selectizeInput("gene",
                   "Gene Name:",
                   choices=c("",sel_list))
  })
  
  output$esssel2 = renderUI({
    sel_list = c("")
    if(input$mutdata == "DEPMAP")
    {
      sel_list = colnames(pdx_plot)[!grepl("^MUT[.]",colnames(pdx_plot))]
    }
    else if(input$mutdata == "TCGADEPMAP(exp.only)")
    {
      sel_list = colnames(pdx_plot)[!grepl("^MUT[.]",colnames(pdx_plot))]
    }
    else if(input$mutdata == "PDXEDEPMAP(exp.only)")
    {
      sel_list = colnames(pdx_plot)[!grepl("^MUT[.]",colnames(pdx_plot))]
    }
    else if(input$mutdata == "TCGADEPMAP(integrated)")
    {
      #sel_list = colnames(tcga_multi_plot)[!grepl("^MUT[.]",colnames(tcga_multi_plot))&!grepl("^DEL[.]",colnames(tcga_multi_plot))]
      sel_list = colnames(tcga_multi_plot)
    }
    sel_list = setdiff(sel_list,"lineage")
    selectizeInput("gene2",
                   "Gene Name:",
                   choices=c("",sel_list))
  })
  
  output$mutsel = renderUI({
    if(input$muttype=="mut"){
      sel_list = c("")
      if(input$mutdata == "DEPMAP")
      {
        sel_list = colnames(depmap_plot)[grepl("^MUT[.]",colnames(depmap_plot))]
      }
      else if(input$mutdata == "TCGADEPMAP(exp.only)")
      {
        #sel_list = colnames(tcga_plot)[grepl("^MUT[.]",colnames(tcga_plot))]
        sel_list = colnames(tcga_mut_data)[grepl("^MUT[.]",colnames(tcga_mut_data))]
      }
      else if(input$mutdata == "PDXEDEPMAP(exp.only)")
      {
        sel_list = colnames(pdx_plot)[grepl("^MUT[.]",colnames(pdx_plot))]
      }
      else if(input$mutdata == "TCGADEPMAP(integrated)")
      {
        #sel_list = colnames(tcga_multi_plot)[grepl("^MUT[.]",colnames(tcga_multi_plot))]
        sel_list = colnames(tcga_mut_data)[grepl("^MUT[.]",colnames(tcga_mut_data))]
      }

      selectizeInput("mutgene",
                     "Mutation:",
                     choices=c("",sel_list))
    }
    else if(input$muttype=="del")
    {
      sel_list = c("")
      if(input$mutdata == "DEPMAP")
      {
        sel_list = colnames(depmap_plot)[grepl("^DEL[.]",colnames(depmap_plot))]
        selectizeInput("mutgene",
                       "Deletion:",
                       choices=c("",sel_list))
      }
      else if(input$mutdata == "TCGADEPMAP(exp.only)")
      {
        sel_list = colnames(tcga_mut_data)[grepl("^DEL[.]",colnames(tcga_mut_data))]
        selectizeInput("mutgene",
                       "Deletion:",
                       choices=c("",sel_list))
      }
      else if(input$mutdata == "PDXEDEPMAP(exp.only)")
      {
        selectizeInput("mutgene",
                       "Deletion data is not avaiable for PDXE",
                       choices=c(""))
      }
      else if(input$mutdata == "TCGADEPMAP(integrated)")
      {
        sel_list = colnames(tcga_mut_data)[grepl("^DEL[.]",colnames(tcga_mut_data))]
        selectizeInput("mutgene",
                       "Deletion:",
                       choices=c("",sel_list))
      }
    }
    else if(input$muttype=="amp")
    {
      sel_list = c("")
      if(input$mutdata == "DEPMAP")
      {
        sel_list = colnames(depmap_plot)[grepl("^AMP[.]",colnames(depmap_plot))]
        selectizeInput("mutgene",
                       "Amplification:",
                       choices=c("",sel_list))
      }
      else if(input$mutdata == "TCGADEPMAP(exp.only)")
      {
        sel_list = colnames(tcga_mut_data)[grepl("^AMP[.]",colnames(tcga_mut_data))]
        selectizeInput("mutgene",
                       "Amplification:",
                       choices=c("",sel_list))
      }
      else if(input$mutdata == "PDXEDEPMAP(exp.only)")
      {
        selectizeInput("mutgene",
                       "Amplification data is not avaiable for PDXE",
                       choices=c(""))
      }
      else if(input$mutdata == "TCGADEPMAP(integrated)")
      {
        sel_list = colnames(tcga_mut_data)[grepl("^AMP[.]",colnames(tcga_mut_data))]
        selectizeInput("mutgene",
                       "Amplification:",
                       choices=c("",sel_list))
      }
    }
    
  })
  
  observeEvent(input$gene,{
    output$distPlot = renderPlot({
      if(input$gene != "")
      {
        x_gene = input$gene
        ggboxplot(xdata(),x="lineage",y=x_gene) + rotate_x_text() + 
          xlab("Lineage") + ylab(paste0(x_gene," essentiality")) + theme(aspect.ratio = 0.5)
      }
    })
    if(input$gene != "")
    {
      output$down = renderUI({
        downloadButton("downlin","Download plot data")
      })
    }
  })
  
  output$downlin <- downloadHandler(
    filename = function() {
      paste('data.csv', sep='')
  },
  content = function(con) {
     write.csv(xdata(), con)
    }
  )
  
  xmutdata = reactive({
    x_data = c("")
    if(input$mutgene != "")
    {
      x_gene = input$mutgene
      if(input$mutdata == "DEPMAP" & input$mutgene%in%colnames(depmap_plot))
      {
        x_data = depmap_plot[,c(x_gene,input$gene2)]
      }
      else if(input$mutdata == "TCGADEPMAP(exp.only)" & input$mutgene%in%colnames(tcga_mut_data))
      {
        x_data = data.frame(tcga_plot,tcga_mut_data[rownames(tcga_plot),])[,c(x_gene,input$gene2)]
      }
      else if(input$mutdata == "PDXEDEPMAP(exp.only)" & input$mutgene%in%colnames(pdx_plot))
      {
        x_data = pdx_plot[,c(x_gene,input$gene2)]
      }
      else if(input$mutdata == "TCGADEPMAP(integrated)" & input$mutgene%in%colnames(tcga_mut_data))
      {
        x_data = data.frame(tcga_multi_plot,tcga_mut_data[rownames(tcga_multi_plot),])[,c(x_gene,input$gene2)]
      }
      if(length(x_data)>1)
      {
        key_id = strsplit(x_gene,"[.]")[[1]][2]
        key_s = "MUT"
        if(input$muttype=="del"){
          key_s = "DEL"
        }
        else if(input$muttype=="amp"){
          key_s = "AMP"
        }
        #key_id = gsub(paste0(key_s,"."),"",x_gene)
        idx0 = which(x_data[,x_gene]==0)
        idx1 = which(x_data[,x_gene]==1)
        x_data[idx0,x_gene] = paste0(key_id,"_WT")
        x_data[idx1,x_gene] = paste0(key_id,"_",key_s)
      }
      x_data
    }
  })
  
  observeEvent(input$mutgene,{
    output$mutPlot = renderPlot({
      if(is.null(input$mutgene)){
        return()
      }
      if(input$mutgene != "" & length(xmutdata())>1)
      {
        x_gene = input$mutgene
        ggboxplot(xmutdata(),x=x_gene,y=input$gene2,outlier.shape = NA) + xlab("") + ylab(paste0(input$gene2," essentiality")) + 
          stat_compare_means() + rotate_x_text(45) + theme(aspect.ratio = 1.5)
      }
    })
    if(input$mutgene != "")
    {
      output$down1 = renderUI({
        downloadButton("downmut","Download plot data")
      })
    }
  })
  
  output$downmut <- downloadHandler(
    filename = function() {
      paste('data.csv', sep='')
    },
    content = function(con) {
      write.csv(xmutdata(), con)
    }
  )
  observeEvent(input$mutdata,{
    removeUI("#downmut")
  })
}

shinyApp(ui = ui, server = server)
