
library(shiny)

source('helpers.R')

gw_autoph_response_data <- load_kinetic_response_data()
gw_autoph_competence_data <- load_autophagy_competence_data()

autophagy_competence_view <- list()
autophagy_competence_view[['single']] <-
  bslib::card(
    full_screen = TRUE,
    #bslib::card_header("Autophagy competence"),
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
    shiny::htmlOutput("gene_info"),
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

      # sidebar = shiny::selectizeInput(
      #     inputId   = "gene_id_kinetic",
      #     label     = "Gene",
      #     choices   = unique(gw_autoph_competence_data$gene_ids),
      # ),

      sidebar = shiny::selectInput(
        "gene_id_kinetic", "Select ORF/Gene mutant",
        gw_autoph_response_data$gene_ids),
      response_kinetics_view[['single']],
      shiny::uiOutput("selected_gene")
    )
  )

kinetics_sidebar_multiple <-
  bslib::page_fillable(
    bslib::layout_sidebar(
      sidebar = shiny::selectizeInput(
        inputId = "gene_id_kinetic_multiple",
        label = "Select ORF/Gene mutants (max 10)",
        selected = NULL,
        options = list(
          minItems = 1,
          maxItems = 10),
        choices = gw_autoph_response_data$gene_ids),
      response_kinetics_view[['multiple']]
    )
  )

about_page <-
  bslib::page_fillable(
    bslib::card(
      #class = "bg-dark",
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
          "gene_id_bf", "Gene",
          gw_autoph_competence_data$gene_ids),
        shiny::selectInput(
          "bf_dnn_model", "Deep neural network (DNN) model",
          c("30","22")
        )
      ),
      autophagy_competence_view[['single']]
    )
  )

ui <- bslib::page_navbar(
  #footer = list(shinyjs::useShinyjs(), shinyauthr::loginUI("login")),
  theme = bslib::bs_theme(
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
    "input-border-color" = "#EA80FC"
  ),
  #theme = bslib::bs_theme(
  #  bootswatch = "united",
  #  danger = "#343a40",
  #),
  title = "Genome-Wide Autophagy Dynamics",
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
    bslib::nav_panel("Single gene analysis", bfactor_sidebar)
    #bslib::nav_panel("Multiple gene analysis", bfactor_sidebar)
  )

  #bslib::nav_item(tags$a("About", href = "https://posit.co")),

)

# Enable thematic
#thematic::thematic_shiny(font = "auto")

# Change ggplot2's default "gray" theme
#ggplot2::theme_set(ggplot2::theme_bw(base_size = 18))

# New server logic (removes the `+ theme_bw()` part)
server <- function(input, output, session) {

  # #login status and info will be managed by shinyauthr module and stores here
  # credentials <- shinyauthr::loginServer(
  #     id = "login",
  #     data = user_base,
  #     user_col = user,
  #     pwd_col = password,
  #     sodium_hashed = T,
  #     log_out = reactive(logout_init())
  # )
  #
  # # call the logout module with reactive trigger to hide/show
  # logout_init <- shinyauthr::logoutServer(
  #     id = "logout",
  #     active = reactive(credentials()$user_auth)
  # )
  #
  # # this opens or closes the sidebar on login/logout
  # observe({
  #     if(credentials()$user_auth) {
  #         shinyjs::removeClass(selector = "body", class = "sidebar")
  #     } else {
  #         shinyjs::addClass(selector = "body", class = "sidebar")
  #     }
  # })

  output$autoph_response_kinetics <- shiny::renderPlot({
    #req(credentials()$user_auth)
    plot_response_kinetics(
      response_data = gw_autoph_response_data$per_ko[[input$gene_id_kinetic]])
  })

  output$autoph_response_kinetics_multiple <- shiny::renderPlot({
    #req(credentials()$user_auth)
    plot_response_kinetics_multiple(
      response_data = gw_autoph_response_data,
      gene_ids = input$gene_id_kinetic_multiple)
  })

  output$autoph_competence <- shiny::renderPlot({
    #req(credentials()$user_auth)
    plot_autophagy_competence(
      competence_data = gw_autoph_competence_data$per_ko[[input$gene_id_bf]],
      dnn_model = input$bf_dnn_model)
  })

  output$selected_gene <- shiny::renderText({
    paste("You have selected", input$gene_id_kinetic)
  })

  output$gene_info <- shiny::renderUI({
    paste("<div><a href='https://www.vg.no'>VG</a></div>")
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
