
library(shiny)

source('helpers.R')

gw_autoph_response_data <- load_kinetic_response_data()
gw_autoph_competence_data <- load_autophagy_competence_data()

for(elem in c('ds_parms','ds_parms_comb')){
  if(!is.null(gw_autoph_response_data[[elem]])){
    gw_autoph_response_data[[elem]]$Parameter <-
      gsub("A_","Autophagy ",gw_autoph_response_data[[elem]]$Parameter)
    gw_autoph_response_data[[elem]]$Parameter <-
      ifelse(!grepl("Param",gw_autoph_response_data[[elem]]$Parameter),
             gsub("slope","Slope (tangent) ",gw_autoph_response_data[[elem]]$Parameter),
             gsub("slope","Slope (sigmoid) ",gw_autoph_response_data[[elem]]$Parameter))
    gw_autoph_response_data[[elem]]$Parameter <- gsub("Param","",gw_autoph_response_data[[elem]]$Parameter)
    gw_autoph_response_data[[elem]]$Parameter <- gsub("1","-N",gw_autoph_response_data[[elem]]$Parameter)
    gw_autoph_response_data[[elem]]$Parameter <- gsub("2","+N ",gw_autoph_response_data[[elem]]$Parameter)
    gw_autoph_response_data[[elem]]$Parameter <- gsub("_"," ",gw_autoph_response_data[[elem]]$Parameter)
    gw_autoph_response_data[[elem]]$Parameter <- gsub("starvation","-N",gw_autoph_response_data[[elem]]$Parameter)
    gw_autoph_response_data[[elem]]$Parameter <- gsub("replenishment","+N",gw_autoph_response_data[[elem]]$Parameter)
    gw_autoph_response_data[[elem]]$Parameter <- stringr::str_trim(gw_autoph_response_data[[elem]]$Parameter)

  }
}

autophagy_competence_view <- list()
autophagy_competence_view[['single']] <-
  bslib::card(
    full_screen = TRUE,
    #bslib::card_header("Autophagy competence"),
    shiny::htmlOutput("gene_info_bf"),
    shiny::plotOutput("autoph_competence")
  )
autophagy_competence_view[['multiple']] <-
  bslib::card(
    full_screen = TRUE,
    #bslib::card_header("Autophagy competence"),
    shiny::plotOutput("autoph_competence_multiple")
  )

response_kinetics_view <- list()
response_kinetics_view[['single']] <-
  bslib::card(
    full_screen = TRUE,
    bslib::card_header(
      shiny::textOutput("selected_gene")),
    shiny::htmlOutput("gene_info_kinetic"),
    #bslib::card_header("Autophagy response kinetics"),
    shiny::plotOutput("autoph_response_kinetics")
  )
response_kinetics_view[['multiple']] <-
  bslib::card(
    full_screen = TRUE,
    #bslib::card_header("Autophagy response kinetics"),
    shiny::plotOutput("autoph_response_kinetics_multiple")
  )

kinetics_sidebar <-
  bslib::page_fillable(
    bslib::layout_sidebar(
      sidebar = shiny::selectInput(
        "gene_id_kinetic", "Select ORF/Gene mutant",
        gw_autoph_response_data$gene_info_kinetic$orf_gene_id, selectize = T),
      response_kinetics_view[['single']],
      shiny::uiOutput("selected_gene")
    )
  )

kinetics_sidebar_multiple <-
  bslib::page_fillable(
    bslib::layout_sidebar(
      sidebar = list(
        shiny::selectizeInput(
        inputId = "gene_id_kinetic_multiple",
        label = "Select ORF/Gene mutants (max 10)",
        selected = "RTG2 / YGL252C",
        options = list(
          minItems = 1,
          maxItems = 10),
        #selectize = T,
        choices = gw_autoph_response_data$gene_info_kinetic$orf_gene_id),
        shiny::selectInput(
          "x_var_kin","X-axis variable",
          c("Perturbation -N",
            "Perturbation +N",
            "Perturbation overall",
            "Slope (sigmoid) -N",
            "Slope (sigmoid) +N",
            "Slope (tangent) -N",
            "Slope (tangent) +N",
            "Autophagy start",
            "Autophagy max",
            "Autophagy final",
            "T50 -N",
            "T50 +N",
            "T lag -N",
            "T lag +N",
            "T final -N",
            "T final +N",
            "Dynamic range -N",
            "Dynamic range +N"),
          selected = "T50 +N"
        ),
        shiny::selectInput(
          "y_var_kin","Y-axis variable",
          c("Perturbation -N",
            "Perturbation +N",
            "Perturbation overall",
            "Slope (sigmoid) -N",
            "Slope (sigmoid) +N",
            "Slope (tangent) -N",
            "Slope (tangent) +N",
            "Autophagy start",
            "Autophagy max",
            "Autophagy final",
            "T50 -N",
            "T50 +N",
            "T lag -N",
            "T lag +N",
            "T final -N",
            "T final +N",
            "Dynamic range -N",
            "Dynamic range +N"),
          selected = "T50 -N"
        ),
        shiny::checkboxInput(
          "contour", "Add contour plot", value = F)),
      response_kinetics_view[['multiple']]
    )
  )

about_page <-
  bslib::page_fillable(
    #home_ui
    bslib::card(
      full_screen = F,
      #bslib::card_header("About", class = "bg-primary text-white"),
      bslib::card_body(
        "This is a web portal for exploring autophagy dynamics in yeast."
      )
    )
  )

downloads_page <-
  bslib::page_fillable(
    bslib::card(
      #class = "bg-dark",
      full_screen = F,
      #bslib::card_header("About", class = "bg-primary text-white"),
      bslib::card_body(
        "Links to file downloads"
      )
    )
  )

bfactor_sidebar <-
  bslib::page_fillable(
    bslib::layout_sidebar(
      sidebar = list(
        shiny::selectInput(
          "gene_id_bf", "Select ORF/Gene mutant",
          gw_autoph_competence_data$gene_info_bf$orf_gene_id),
        shiny::selectInput(
          "bf_dnn_model", "Deep neural network (DNN) model",
          c("30","22")
        )
      ),
      autophagy_competence_view[['single']]
    )
  )


bfactor_sidebar_multiple <-
  bslib::page_fillable(
    bslib::layout_sidebar(
      sidebar = list(
        shiny::selectizeInput(
          inputId = "gene_id_bf_multiple",
          label = "Select ORF/Gene mutants (max 10)",
          selected = "RTG2 / YGL252C",
          choices =
            gw_autoph_competence_data$gene_info_bf$orf_gene_id,
          options = list(
            minItems = 1,
            maxItems = 10)),
        #shiny::selectInput(
        #  "gene_id_bf_multiple", "Gene",
        #  gw_autoph_competence_data$gene_info_bf$orf_gene_id),
        shiny::selectInput(
          "x_var","X-axis variable",
          c("Overall autophagy",
            "Autophagosome formation",
            "Autophagosome clearance"),
          selected = "Autophagosome formation"
        ),
        shiny::selectInput(
          "y_var","Y-axis variable",
          c("Overall autophagy",
            "Autophagosome formation",
            "Autophagosome clearance"),
          selected = "Autophagosome clearance"
        ),
        shiny::selectInput(
          "bf_dnn_model_multiple", "Deep neural network (DNN) model",
          c("30","22")
        )
      ),
      autophagy_competence_view[['multiple']]
    )
  )

page_bs_theme <-
  bslib::bs_theme(
    # Controls the default grayscale palette
    bg = "#ffffff",
    fg = "black",
    # Controls the accent (e.g., hyperlink, button, etc) colors
    primary = "black", secondary = "#48DAC6",
    base_font = c("Grandstander", "sans-serif"),
    code_font = c("Courier", "monospace"),
    bootswatch = "bootstrap",
    heading_font = "'Helvetica Neue', Helvetica, sans-serif",
    # Can also add lower-level customization
    "input-border-color" = "#EA80FC")
  # bslib::bs_add_rules(
  #   rules = "navbar navbar-default {
  #                       height: 400px;
  #                       background-color: pink !important;
  #               }"
  # )

ui <- bslib::page_navbar(
  #footer = list(shinyjs::useShinyjs(), shinyauthr::loginUI("login")),
  theme = page_bs_theme,
  #theme = bslib::bs_theme(
  #  bootswatch = "united",
  #  danger = "#343a40",
  #),
  title = "Genome-Wide Autophagy Dynamics in Yeast",
  #subtitle = "A web portal for exploring autophagy dynamics in yeast",
  bslib::nav_spacer(),
  bslib::nav_panel("Home", about_page),
  bslib::nav_panel("Downloads", downloads_page),
  bslib::nav_menu(
    "Autophagy response kinetics",
    bslib::nav_panel("Single gene perspective", kinetics_sidebar),
    bslib::nav_panel("Multiple genes perspective", kinetics_sidebar_multiple)
  ),
  bslib::nav_menu(
    "Autophagy competence",
    bslib::nav_panel("Single gene perspective", bfactor_sidebar),
    bslib::nav_panel("Multiple gene perspective", bfactor_sidebar_multiple)
  )

  #bslib::nav_item(tags$a("About", href = "https://posit.co")),

)

# Enable thematic
#thematic::thematic_shiny(font = "auto")

# Change ggplot2's default "gray" theme
#ggplot2::theme_set(ggplot2::theme_bw(base_size = 18))

# New server logic (removes the `+ theme_bw()` part)
server <- function(input, output, session) {

  output$autoph_response_kinetics <- shiny::renderPlot({
    plot_response_kinetics(
      response_data = gw_autoph_response_data$per_ko[[input$gene_id_kinetic]])
  })

  output$autoph_response_kinetics_multiple <- shiny::renderPlot({
    plot_response_kinetics_multi(
      response_data = gw_autoph_response_data,
      user_x = input$x_var_kin,
      user_y = input$y_var_kin,
      show_library_type_contour = input$contour,
      primary_identifiers = input$gene_id_kinetic_multiple)
  })

  output$autoph_competence_multiple <- shiny::renderPlot({
    plot_autophagy_competence_multi(
      competence_data = gw_autoph_competence_data,
      dnn_model = input$bf_dnn_model_multiple,
      user_x = input$x_var,
      user_y = input$y_var,
      primary_identifiers = input$gene_id_bf_multiple)
  })

  output$autoph_competence <- shiny::renderPlot({
    plot_autophagy_competence(
      competence_data = gw_autoph_competence_data$per_ko[[input$gene_id_bf]],
      dnn_model = input$bf_dnn_model)
  })

  output$selected_gene <- shiny::renderText({
    paste(input$gene_id_kinetic)
  })

  output$gene_info_kinetic <- shiny::renderUI({
    ginf <- show_gene_info(
      primary_id = input$gene_id_kinetic,
      gene_info = gw_autoph_response_data$gene_info_kinetic
    )
    shiny::HTML("<div><ul><li>Genename: ",ginf[['sgd_link']],"</li>",
                "<li>Description: ",ginf[['description']],"</li>",
                "<li>Human orthologs: ",ginf[['human_orthologs']],"</li>",
                "</ul></div><br>")
  })

  output$gene_info_bf <- shiny::renderUI({
    ginf <- show_gene_info(
      primary_id = input$gene_id_bf,
      gene_info = gw_autoph_competence_data$gene_info_bf
    )
    shiny::HTML("<div><ul><li>Genename: ",ginf[['sgd_link']],"</li>",
                "<li>Description: ",ginf[['description']],"</li>",
                "<li>Human orthologs: ",ginf[['human_orthologs']],"</li>",
                "</ul></div><br>")
  })

  # observe({
  #     updateSelectizeInput(
  #         session,
  #         inputId = 'gene_id_bf',
  #         label   = 'Gene',
  #         selected = head(unique(gw_autoph_competence_data$gene_ids), 1),
  #         choices = unique(gw_autoph_competence_data$gene_ids)
  #     )
  # })
}

# Run the application
shiny::shinyApp(ui = ui, server = server)
