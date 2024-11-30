
library(shiny)

source('helpers.R')

gw_autoph_response_data <- load_kinetic_response_data()
gw_autoph_competence_data <- load_autophagy_competence_data()
kinetic_response_params <- get_kinetic_response_params()

navbar_js <- "@media (max-width: 800px) {
    .navbar-header {
        float: none;
    }
    .navbar-left,.navbar-right {
        float: none !important;
    }
    .navbar-toggle {
        display: block;
    }
    .navbar-collapse {
        border-top: 1px solid transparent;
        box-shadow: inset 0 1px 0 rgba(255,255,255,0.1);
    }
    .navbar-fixed-top {
        top: 0;
        border-width: 0 0 1px;
    }
    .navbar-collapse.collapse {
        display: none!important;
    }
    .navbar-nav {
        float: none!important;
        margin-top: 7.5px;
    }
    .navbar-nav>li {
        float: none;
    }
    .navbar-nav>li>a {
        padding-top: 10px;
        padding-bottom: 10px;
    }
    .collapse.in{
        display:block !important;
    }
}"

autophagy_competence_view <- list()
autophagy_competence_view[['single']] <-
  bslib::page_fillable(
    bslib::card(
      full_screen = TRUE,
      #bslib::card_header("Autophagy competence"),
      shiny::htmlOutput("gene_info_bf"),
      shiny::plotOutput("autoph_competence")
    ),
    shiny::includeHTML("data/section_content/citation_footnote.md")
  )

autophagy_competence_view[['multiple']] <-
  bslib::card(
    full_screen = TRUE,
    shiny::plotOutput("autoph_competence_multiple")
  )

response_kinetics_view <- list()
response_kinetics_view[['single']] <-
  bslib::page_fillable(
    bslib::card(
      full_screen = TRUE,
      bslib::card_header(
        shiny::textOutput("selected_gene")),
      shiny::htmlOutput("gene_info_kinetic"),
      shiny::plotOutput("autoph_response_kinetics")
    ),
    shiny::includeHTML("data/section_content/citation_footnote.md")
  )
response_kinetics_view[['multiple']] <-
  bslib::card(
    full_screen = TRUE,
    shiny::plotOutput("autoph_response_kinetics_multiple")
  )

kinetics_sidebar <-
  bslib::page_fillable(
    bslib::layout_sidebar(
      sidebar = shiny::selectInput(
        "gene_id_kinetic", "Select gene/ORF mutant",
        gw_autoph_response_data$gene_info_kinetic$orf_gene_id,
        width = "100%"),
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
        label = "Highlight gene/ORF mutants (max 10)",
        selected = "AAC1 / YMR056C",
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
  "<h5><b><p style=text-align:justify;'>",
  "<a href='https://www.biorxiv.org/content/10.1101/2024.04.06.588104v1' ",
  "target='_blank'>Genome-wide profiling of the hierarchical control of autophagy ",
  "dynamics using deep learning (bioRxiv, 2024)</a></p></b></h5>",
  "<h6><p style='text-align:justify;'>",
  htmltools::includeText("data/section_content/authors.md"),
  "</p></h6>",
  "<p style='text-align:justify;'><h6>Correspondence to <i>nathac@uio.no </i>",
  "or <i>j.m.enserink@ibv.uio.no</i></h6></p>")

synopsis <- paste0(
  "<h5><b><p style='color:#593196;text-align:justify;'>",
  htmltools::includeText("data/section_content/synopsis_I.md"),
  "</p></b></h5>",
  "<h6><p style='text-align:justify;'>",
  htmltools::includeText("data/section_content/synopsis_II.md"),
  "</p></h6>")

what_is_autodry <- paste0(
  "<h5><b><p style='color:#593196;text-align:justify;'>",
  htmltools::includeText("data/section_content/what_is_autodry_I.md"),
  "</p></b></h5>",
  "<h6><p style='text-align:justify;'>",
  htmltools::includeText("data/section_content/what_is_autodry_II.md"),
  "</p></h6>")

about_study_text <- paste0(
  "<p style='text-align:justify;'>",
  htmltools::includeText("data/section_content/about_the_study_I.md"),
  "</p>",
  "AutoDRY provides two key outputs that provide different ",
  "perspectives on autophagy regulation:<br>",
  "<p style='text-align:justify;'>",
  htmltools::includeText("data/section_content/about_the_study_II.md"),
  "</p>")

disclaimer_text <- paste0(
  "<p style='text-align:justify;'>",
  htmltools::includeText("data/section_content/disclaimer.md"),
  "</p>")

autodry_footer <-
  bslib::card_footer(
    shiny::markdown(
      paste0(
        "<br>",
        "<div align='center'>",
        "<a href='https://www.uio.no' target='_blank'><img src='uio.png' alt='uio' style='width:20%;height:90%;'></a>",
        "<a href='https://ous-research.no/institute' target='_blank'><img src='ous.png' alt='ous' style='width:20%'></a>",
        "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",
        "<a href='https://www.med.uio.no/cancell/english/' target='_blank'><img src='cancell.png' style='width:4%;height:48%;'></a>",
        "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",
        "<a href='https://arctic-autophagy.no/' target='_blank'><img src='arctic.png' style='width:10%;height:45%'></a>",
        "</div>")
    )
  )

about_page <-
  bslib::page_fillable(
    tags$style(HTML(".card {border-radius: 0.9rem;}")),
    tags$style(HTML(".card2 {border-radius: 0.1rem;}")),
    bslib::card(
      class = "card2",
      full_screen = F,
      fillable = F,
      fill = F,
      bslib::card_body(
        #min_height = 300,
        bslib::layout_column_wrap(
          width = NULL,
          fill = F,
          style = bslib::css(
            grid_template_columns = "22fr 1fr 22fr 1fr 22fr"),
          bslib::card(
            bslib::card_header(class = "bg-dark","  Background"),
            shiny::markdown(synopsis)),
          shiny::markdown(""),
          bslib::card(
            bslib::card_header(class = "bg-dark","  What is AutoDRY?"),
            shiny::markdown(what_is_autodry)),
          shiny::markdown(""),
          bslib::card(
            bslib::card_header(class = "bg-dark","  Citation"),
            shiny::markdown(paper_info)
          )
        ),
        shiny::markdown(
          paste0(
            "<br><br><div align='center'>",
            "<img src='MNS_Figure_1A_v2.png' alt='Autophagy_Dynamics_Study_Overview'",
            " width='90%' height='100%' align='center'/></div><br>"
          )
        ),
        bslib::layout_column_wrap(
          width = NULL,
          fill = F,
          bslib::card(
            bslib::card_header(class = "bg-dark","  About the study"),
            shiny::markdown(about_study_text)
          )
        )
      ),
      autodry_footer
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

disclaimer_page <-
  bslib::page_fillable(
    tags$style(HTML(".card {border-radius: 0.9rem;}")),
    tags$style(HTML(".card2 {border-radius: 0rem}")),
    bslib::card(
      class = "card2",
      full_screen = F,
      fillable = F,
      fill = F,
      bslib::card_body(
        bslib::layout_column_wrap(
          width = NULL,
          fill = T,
          bslib::card(
            bslib::card_header(class = "bg-dark","  DISCLAIMER"),
            bslib::card_body(
              shiny::markdown(disclaimer_text)
            )
          )
        )
      ),
      autodry_footer
    )
  )


bfactor_sidebar <-
  bslib::page_fillable(
    bslib::layout_sidebar(
      sidebar = list(
        shiny::selectInput(
          "gene_id_bf", "Select gene/ORF mutant",
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
          label = "Highlight gene/ORF mutants (max 10)",
          selected = "AAC1 / YMR056C",
          choices =
            gw_autoph_competence_data$gene_info_bf$orf_gene_id,
          options = list(
            minItems = 1,
            maxItems = 10)),
        #shiny::selectInput(
        #  "gene_id_bf_multiple", "Gene",
        #  gw_autoph_competence_data$gene_info_bf$orf_gene_id),
        shiny::selectInput(
          "bf_x_var",
          "X-axis variable",
          c("Overall autophagy",
            "Autophagosome formation",
            "Autophagosome clearance"),
          selected = "Autophagosome formation"
        ),
        shiny::selectInput(
          "bf_y_var",
          "Y-axis variable",
          c("Overall autophagy",
            "Autophagosome formation",
            "Autophagosome clearance"),
          selected = "Autophagosome clearance"
        ),
        shiny::selectInput(
          "bf_dnn_model_multiple",
          "Deep neural network (DNN) model",
          c("30","22")
        ),
        shiny::checkboxInput(
         "bf_library_adjustment", "Perform library correction", value = F),
        shiny::checkboxInput(
          "bf_contour", "Add contour plot", value = F)
      ),
      autophagy_competence_view[['multiple']]
    )
  )


ui <- bslib::page_navbar(
  tags$head(
    shiny::includeHTML("google_analytics.html"),
    tags$link(rel="shortcut icon", href="favicon-ous.svg")),
  # tags$style(HTML(
  #  '
  #      .navbar-brand {font-size:1.2em;}
  #      .nav.navbar-nav {font-size:1.2em;}
  #     '
  # )),
  #footer = list(shinyjs::useShinyjs(), shinyauthr::loginUI("login")),
  #theme = page_bs_theme,
  fillable_mobile = TRUE,
  theme = bslib::bs_theme(
    bootswatch = "pulse") |>
    bslib::bs_add_rules(
      list(
        ".navbar {padding-left: 10px; padding-right: 20px;}",
        ".nav.navbar-nav {font-size:1.2em;}",
        ".navbar.navbar-default {
          background-color: $primary !important;
        }"
      )
    ),
  window_title = "AutoDRY",
  title = htmltools::span(
    "AutoDRY: Autophagy Dynamics Repository Yeast",
    style="font-size:1.6rem; padding-left:10px;"),
  bslib::nav_spacer(),
  bslib::nav_panel("Home", about_page),
  bslib::nav_menu(
    "Data Explorer",
    bslib::nav_panel("Response kinetics - individual", kinetics_sidebar),
    bslib::nav_panel("Response kinetics - global", kinetics_sidebar_multiple),
    bslib::nav_panel("Autophagy competence - individual", bfactor_sidebar),
    bslib::nav_panel("Autophagy competence - global", bfactor_sidebar_multiple)
  ),
  # bslib::nav_menu(
  #   "Autophagy response kinetics",
  #   bslib::nav_panel("Single gene perspective", kinetics_sidebar),
  #   bslib::nav_panel("Multiple genes perspective", kinetics_sidebar_multiple)
  # ),
  # bslib::nav_menu(
  #   "Autophagy competence",
  #   bslib::nav_panel("Single gene perspective", bfactor_sidebar),
  #   bslib::nav_panel("Multiple gene perspective", bfactor_sidebar_multiple)
  # ),
  bslib::nav_panel("Downloads", downloads_page),
  bslib::nav_panel(htmltools::span(
    "DISCLAIMER", style="padding-right:10px;"), disclaimer_page)
  #bslib::nav_panel("DISCLAIMER", disclaimer_page)
)

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
      user_x = input$bf_x_var,
      user_y = input$bf_y_var,
      show_library_type_contour = input$bf_contour,
      library_adjustment = input$bf_library_adjustment,
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
    shiny::HTML(
      "<div><ul><li>Genename: ",ginf[['sgd_link']],"</li>",
      "<li>Description: ",ginf[['description']],"</li>",
      "<li>Human orthologs: ",ginf[['human_orthologs']],"</li>",
      "<li>Autophagy perturbation response profile: ",
      ginf[['response_profile']],"</li>",
      "</ul></div>")
  })

  output$gene_info_bf <- shiny::renderUI({
    ginf <- show_gene_info(
      primary_id = input$gene_id_bf,
      gene_info = gw_autoph_competence_data$gene_info_bf
    )
    shiny::HTML(
      "<div><ul><li>Genename: ",ginf[['sgd_link']],"</li>",
      "<li>Description: ",ginf[['description']],"</li>",
      "<li>Human orthologs: ",ginf[['human_orthologs']],"</li>",
      "</ul></div>")
  })

}

# Run the application
shiny::shinyApp(ui = ui, server = server)
