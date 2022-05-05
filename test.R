library(shiny)
library(tools)

runApp(
  list(
    ui = tabPanel("File Upload",
                  h4("File Upload"),
                  sidebarLayout(
                    sidebarPanel(
                      fileInput('file1', 'Choose CSV file to upload',
                                accept = c(
                                  'text/csv',
                                  'text/comma-separated-values',
                                  'text/tab-separated-values',
                                  'text/plain',
                                  'csv',
                                  'tsv'
                                )
                      )
                    ),
                    mainPanel(
                      tableOutput('upload')
                    )
                  )
                  
    ),
    server = function(input, output){
      output$upload <- renderTable({
        
        #assign uploaded file to a variable
        File <- input$file1        
        
        #catches null exception
        if (is.null(File))
          return(NULL)
        
        validate(
          need(file_ext(File$name) %in% c(
            'text/csv',
            'text/comma-separated-values',
            'text/tab-separated-values',
            'text/plain',
            'csv',
            'tsv'
          ), "Wrong File Format try again!"))
        
        read.csv(File$datapath)
      })
    }
  ))
