#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

source('helpers.R')

gw_autoph_response_data <- load_kinetic_response_data()
gw_autoph_competence_data <- load_autophagy_competence_data()

gw_plots <- list(
    bslib::card(
        full_screen = TRUE,
        bslib::card_header("Autophagy competence"),
        shiny::plotOutput("response_competence")
    ),
    bslib::card(
        full_screen = TRUE,
        bslib::card_header("Response kinetics"),
        shiny::plotOutput("response_kinetics")
    )
)

kinetics_sidebar <-
    bslib::page_fillable(
        bslib::layout_sidebar(
            sidebar = shiny::selectInput(
                "gene_id_kinetic", "Select gene",
                gw_autoph_response_data$gene_ids),
            gw_plots[[2]]
        )
    )

bfactor_sidebar <-
    bslib::page_fillable(
        bslib::layout_sidebar(
            sidebar = shiny::selectInput(
                "gene_id_bf", "Select gene",
                gw_autoph_competence_data$gene_ids),
            gw_plots[[1]]
        )
    )
ui <- bslib::page_navbar(
    title = "Genome-Wide Autophagy Dynamics",
    #subtitle = "A web portal for exploring autophagy dynamics in yeast",
    bslib::nav_spacer(),
    bslib::nav_panel("Response kinetics", kinetics_sidebar),
    bslib::nav_panel("Autophagy competence", bfactor_sidebar),
    bslib::nav_item(tags$a("About", href = "https://posit.co")),
    theme = bslib::bs_theme(
        bootswatch = "united"
    )
)

# Enable thematic
#thematic::thematic_shiny(font = "auto")

# Change ggplot2's default "gray" theme
#ggplot2::theme_set(ggplot2::theme_bw(base_size = 18))

# New server logic (removes the `+ theme_bw()` part)
server <- function(input, output) {

    output$response_kinetics <- shiny::renderPlot({
        plot_response_kinetics(response_data = gw_autoph_response_data$per_ko[[input$gene_id_kinetic]])
     })

    output$response_competence <- shiny::renderPlot({
        plot_autophagy_competence(competence_data = gw_autoph_competence_data$per_ko[[input$gene_id_bf]])
    })
}

# Run the application
shinyApp(ui = ui, server = server)
