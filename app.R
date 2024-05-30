library(shiny)
options(rsconnect.packrat=TRUE)
source('helpers.R')

# jsonlite problem here
# https://forum.posit.co/t/error-deploying-jsonlite-fails-to-build/186252

server <- function(input, output) {

    grid <- seq(15, 50, by = .5)

    L <- reactive({
        choleski.decompose.SE(
            alpha=input$alpha,
            rho=input$rho,
            grid=grid
        )
    })

    draws <- reactive({ 
        lapply(1:input$n_realisations, function(x) {
            l <- length(grid)
            rnorm(l)
        })
    })

    realisations <- reactive({
        draws.to.gp(mu = input$mu, normal_draws = draws(), L = L(), grid = grid)
    })

    realisations_logit <- reactive({
        realisations()[, y := inverse.logit(y)]
    })

    output$gp_plot <- renderPlot({
        plot.gaussian.process.realisations(
            draws = realisations(),
            mu = input$mu,
            type = "base",
            ylims = c(-5,5)
        )
    })

    output$gp_plot_logit <- renderPlot({
        plot.gaussian.process.realisations(
            draws = realisations_logit(),
            mu = input$mu,
            type = "logit",
            ylims = c(0,1)
        )
    })

    output$invgamma_pdf_plot <- renderPlot({
        plot.inverse.gamma(rate = input$ig_rate, shape = input$ig_shape)
    })
}

# select a random bootstrap theme (for exploration)
# selected_theme <- sample(bslib::bootswatch_themes(), size = 1)
selected_theme <- "lux"

# Define UI for app
ui <- shinyUI(fluidPage(
    theme = bslib::bs_theme(bootswatch = selected_theme),

    # App title ---
    titlePanel("GP visualisation: the Squared Exponential kernel."),

    # Sidebar anel for inputs ---
    sidebarLayout(
        sidebarPanel(
            width = 3,
            conditionalPanel(
                "input$tabs1 == 'gpreal'",
                sliderInput(
                    inputId = "mu",
                    label = "Mean mu",
                    min = -5,
                    max = 5,
                    value = 0
                ),
                sliderInput(
                    inputId = "alpha",
                    label = "Marginal standard deviation alpha",
                    min = .1,
                    max = 2,
                    value = 0.2
                ),
                sliderInput(
                    inputId = "rho",
                    label = "Lengthscale rho",
                    min = 0,
                    max = 15,
                    value = 5,
                ),
                sliderInput(
                    inputId = "n_realisations",
                    label = "Number of gp realisations",
                    min = 1,
                    max = 50,
                    value = 5
                )
            ),
            conditionalPanel(
                "input$tabs1 == 'igprior'",
                sliderInput(
                    inputId = "ig_rate",
                    label = "Rate of inverse gamma prior",
                    min = 0,
                    max = 10,
                    value = 2,
                ),
                sliderInput(
                    inputId = "ig_shape",
                    label = "Shape of inverse gamma prior",
                    min = 0,
                    max = 10,
                    value = 2,
                ),
            )
        ),
        mainPanel(
            width = 8,
            p(" This webpage aims to provide intuitions on how functions sampled from Gaussian Processes are affected by their hyperparameters."),
            p("One of the most common kernels is the squared exponential kernel, which is parametrised by the lengthscale (which I denote by the greek letter $rho$) and the marginal variance (denoted by $alpha$)."),
            p("I consider a slightly more complicated scenario whereby the function I want to model represents a probability of success, and hence has to lie in [0,1]. As such, we take an inverse logit to map the Gaussian Process samples to [0,1]."),
            p("Happy exploration!"
            ),
            tabsetPanel(
                type = "tabs", id = "tabs1",
                tabPanel(
                    value = "gpreal",
                    id = "tab_gp_realisations",
                    title = "GP realisations",
                    plotOutput(outputId = "gp_plot")
                ),
                tabPanel(
                    value = "gpreal",
                    id = "tab_gp_realisations_logit",
                    title = "logit(GP realisations)",
                    plotOutput(outputId = "gp_plot_logit")
                ),
                tabPanel(
                    value = "igprior",
                    id = "tab_ig_prior",
                    title = "Inverse Gamma prior",
                    plotOutput(outputId = "invgamma_pdf_plot")
                ),
            )
        )
    ),
    # main panel for displaying outputs

    p(em(paste("The bootstrap theme used for this webpage is:", selected_theme)))
))

shinyApp(ui=ui, server=server)
