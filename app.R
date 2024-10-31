
library(shiny)

source('helpers.R')

gw_autoph_response_data <- load_kinetic_response_data()
gw_autoph_competence_data <- load_autophagy_competence_data()
kinetic_response_params <- get_kinetic_response_params()

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
        choices = gw_autoph_response_data$gene_info_kinetic_multi$orf_gene_id),
        shiny::selectInput(
          "x_var_kin","X-axis variable",
          kinetic_response_params,
          selected = "T50 +N"
        ),
        shiny::selectInput(
          "y_var_kin","Y-axis variable",
          kinetic_response_params,
          selected = "T50 -N"
        ),
        shiny::checkboxInput(
          "use_perturbation_data", "Use normalized values", value = F),
        shiny::checkboxInput(
          "contour", "Add contour plot", value = F)),
      response_kinetics_view[['multiple']]
    )
  )


paper_info <- paste0(
  "<br><h6><i><a href='https://www.biorxiv.org/content/10.1101/2024.04.06.588104v1' target='_blank'>Genome-wide profiling of the hierarchical control of autophagy dynamics using deep learning</a></i></h6>",
  "Nathalia Chica, Aram N. Andersen, Sara Orellana-Muñoz, Ignacio Garcia, Aurélie Nguéa P, ",
  "Pilar Ayuda-Durán, Linda Håkensbakken, Eline Rødningen, Christopher D. Putnam, ",
  "Manuela Zucknick, Tor Erik Rusten and Jorrit M. Enserink. ",
  "Correspondence: nathac@uio.no or j.m.enserink@ibv.uio.no")

sales_pitch <- paste0(
  "<br><b><p style='color:#593196;'>How is autophagy regulated over time and space, and how can we effectively ",
  "harness the mechanisms that modulate this essential process? While the ",
  "autophagy field has provided detailed insights into core molecular ",
  "regulatory pathways, a lack of systems-level data on autophagy ",
  "control and dynamics has constrained the development of predictive ",
  "models and in vivo manipulation strategies. We contend that answering ",
  "these questions requires a comprehensive understanding of genome-wide ",
  "influences on autophagy dynamics, beyond basic on-off mechanisms.</p></b><br><br>")

about_study_text <- paste0(
  "<h6>About Our Study</h6>",
  "We engineered a yeast library of nearly 6,000 gene deletions, ",
  "each tagged with fluorescent markers that allow us to visualize ",
  "autophagy as it unfolds. By placing these yeast strains under ",
  "conditions of nitrogen starvation for 12 hours followed by ",
  "nitrogen replenishment for 8 hours, we captured hourly snapshots ",
  "of cellular responses, creating a time-resolved, high-content dataset.<br><br>",
  "To make sense of this immense dataset, we applied deep learning to ",
  "predict which cells are actively undergoing autophagy. This allowed us ",
  "to identify subtle differences in autophagy activity and map these ",
  "dynamics onto a “feature space” that captures each gene's unique ",
  "impact on autophagy over time.<br><br>",
  "Our analysis generates two key outputs that provide ",
  "different perspectives on autophagy regulation:
  <ol><li><b>Nutrient Response Kinetics:</b> By tracking the percentage of ",
  "cells with active autophagy, we quantified how each mutant strain responds ",
  "to changes in nutrient availability. Using a double-sigmoidal curve ",
  "fitting approach, we captured the kinetics of autophagy ",
  "activation and deactivation. This analysis reveals the sensitivity and ",
  "timing of each gene’s response, showing how quickly autophagy ramps up ",
  "during nutrient starvation and subsides when nutrients are replenished.</li>",
  "<li><b>Latent Space Analysis of Autophagic Stages:</b> In addition to ",
  "response kinetics, we applied Bayesian analysis to the neural network’s ",
  "latent space features, identifying distinct autophagic profiles across ",
  "mutants. This approach highlights genes with specific disruptions in ",
  "either the formation or clearance of autophagosomes, giving us a window ",
  "into the multi-stage regulation of autophagy. By pinpointing which genes ",
  "impact these distinct stages, we gain deeper insights into how autophagy ",
  "is finely controlled at each step of the process.</li></ol><br>",
  "Through this web portal, you can dive into the autophagy profiles of ",
  "individual genes or analyze groups of genes in the context of ",
  "genome-wide patterns. Both nutrient response kinetics and Bayesian ",
  "analysis of autophagic stages are available, providing insights into ",
  "each gene’s activation timing and specific roles in ",
  "autophagosome formation and clearance.")


about_page <-
  bslib::page_fillable(
    #home_ui
    bslib::card(
      full_screen = F,
      fillable = F,
      fill = F,
      #bslib::card_header("About", class = "bg-primary text-white"),
      bslib::card_body(
        #min_height = 300,
        bslib::layout_column_wrap(
          width = 1/2,
          shiny::markdown(sales_pitch),
          shiny::markdown(paper_info)
        )
      ),
      bslib::card_body(
        shiny::markdown(
          paste0(
            "<div align='center'><img src='MNS_Figure_1A_v2.png' alt='Autophagy_Dynamics_Study_Overview' width='85%' height='95%' align='center'/></div>",
            about_study_text,
            "<hr>",
            "<div align='center'><img src='uio.png' alt='uio' width='20%' height='90%'/>",
            "<img src='ous.png' width='20%'/>",
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src='cancell2.png' width='5%' height='50%'/>",
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src='arctic.png' width='8%' height='50%'/></div>")
        )
      )
      # bslib::card_body(
      #   fillable = F,
      #   #min_height = 150,
      #   max_height = "150px",
      #   height = "150px",
      #   shiny::markdown(
      #     paste0(
      #       "<hr><br>",
      #       "<div align='center'><img src='uio.png' alt='uio' width='20%' height='90%'></div>")
      #   ),
      #   fill = T
      # )
    )
  )

downloads_page <-
  bslib::page_fillable(
    bslib::card(
      #class = "bg-dark",
      full_screen = F,
      #bslib::card_header("About", class = "bg-primary text-white"),
      bslib::card_body(
        shiny::markdown("Links to file downloads will be available here.")
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
  #bslib::bs_theme(preset = "bootstrap")
  bslib::bs_theme(bootswatch = "united")
page_bs_theme <-
  bslib::bs_theme_update(theme = page_bs_theme, bg = "pink", fg = "black", primary="green")

    # # Controls the default grayscale palette
    # bg = "#ffffff",
    # fg = "black",
    # # Controls the accent (e.g., hyperlink, button, etc) colors
    # primary = "black",
    # secondary = "#48DAC6",
    # base_font = c("Grandstander", "sans-serif"),
    # code_font = c("Courier", "monospace"),
    # bootswatch = "bootstrap",
    # heading_font = "'Helvetica Neue', Helvetica, sans-serif") |>
  # bslib::bs_add_rules(
  #   rules = ".navbar .navbar-default .navbar-inverse .navbar-static-top {
  #                       height: 400px;
  #                       color: pink !important;
  #                       background-color: pink !important;
  #               }"
  # )

ui <- bslib::page_navbar(
  #footer = list(shinyjs::useShinyjs(), shinyauthr::loginUI("login")),
  #theme = page_bs_theme,
  theme = bslib::bs_theme(bootswatch = "pulse") |>
    bslib::bs_add_rules(
      rules = "
          .navbar.navbar-default {
                background-color: $primary !important;
          }
          "
    ),
  #theme = bslib::bs_theme(
  #  bootswatch = "united",
  #  danger = "#343a40",
  #),
  title = paste0(
    "AutoDRY: Genome-Wide Autophagy Dynamics Repository Yeast"),
  #subtitle = "A web portal for exploring autophagy dynamics in yeast",
  bslib::nav_spacer(),
  bslib::nav_panel("Home", about_page),
  bslib::nav_menu(
    "Autophagy response kinetics",
    bslib::nav_panel("Single gene perspective", kinetics_sidebar),
    bslib::nav_panel("Multiple genes perspective", kinetics_sidebar_multiple)
  ),
  bslib::nav_menu(
    "Autophagy competence",
    bslib::nav_panel("Single gene perspective", bfactor_sidebar),
    bslib::nav_panel("Multiple gene perspective", bfactor_sidebar_multiple)
  ),
  bslib::nav_panel("Data downloads", downloads_page),

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
      use_perturbation_data = input$use_perturbation_data,
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
