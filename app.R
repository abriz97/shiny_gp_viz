# shiny::runApp("app.R")
library(shiny)
library(ggplot2)
library(data.table)
library(invgamma)
library(MASS)
library(viridis)

# plot.range.alpha <- function(x_range = seq(-9, 9, by = .5), delta = 1) {
#     df_base <- data.table(
#         x = seq(-10, 10, by = 1 / 1000),
#         y = sapply(seq(-10, 10, by = 1 / 1000), inverse.logit)
#     )
#     df <- data.table(
#         x = x_range,
#         y = sapply(x_range, inverse.logit),
#         ymin = sapply(x_range - delta / 2, inverse.logit),
#         ymax = sapply(x_range + delta / 2, inverse.logit)
#     )
#     ggplot(df_base, aes(x = x, y = y)) +
#         geom_line() +
#         geom_point(data = df) +
#         geom_linerange(data = df, aes(ymin = ymin, ymax = ymax, color = y)) +
#         theme_minimal() +
#         labs(y = "inverse logit:\nexp(x)/(1+exp(x)) ")
# }

inverse.logit <- function(x) {
    exp(x) / (1 + exp(x))
}

plot.inverse.gamma <- function(rate = input$ig_rate,
                               shape = input$ig_shape) {
    x_range <- seq(0.05, 50, by = .05)
    df <- data.frame(
        x = x_range,
        y = invgamma::dinvgamma(x = x_range, rate = rate, shape = shape)
    ) 
    # subset to not to small values
    df <- df[df$y >= 1 / 10000, ]
    ggplot(df, aes(x = x, y = y)) +
        geom_line(color = "darkred") +
        theme_minimal() +
        labs(y = "Inverse Gamma pdf")
}

choleski.decompose.SE <- function(alpha, rho, grid){
    kSE <- function(x1, x2) {
        alpha^2 * exp(-(x1 - x2)^2 / (2 * rho^2))
    }
    Sigma <- lapply(grid, kSE, x2 = grid) |> Reduce(f = rbind)
    L <- chol(Sigma + 1e-6 * diag(nrow(Sigma)))
    return(L)
}

draws.to.gp <- function(mu=input$mu, normal_draws=draws, L=L, grid){
    # returns a data table of dimension (n_realisations, length(grid))
    # each row is a single realisation of the Gaussian Process
    gp_realisation <- lapply(normal_draws, function(z){
        out <- z %*% L + mu 
        as.data.table(out)
    }) |> rbindlist()
    names(gp_realisation) <- as.character(grid)
    gp_realisation[, ID := as.factor(1:.N)]
    gp_realisation <- melt.data.table(gp_realisation,
        id.vars = "ID", variable.name = "x", value.name = "y"
    )
    gp_realisation[, x := as.numeric(x)]
    return(gp_realisation)
}

plot.gaussian.process.realisations <- function(draws, mu, type="base", ylims=c(-5,5)){

    # check args
    type <- match.arg(type, c("base", "logit"))
    scale_y <- if(type == "logit") scales::percent else scales::number

    ggplot(draws, aes(x = x, y = y, color = ID, group = ID)) +
        geom_line() +
        geom_label(
            data = data.frame(),
            x = 15.5, y = .95, color = "black", group = NA,
            label = paste0("Mu=", mu)
        ) +
        scale_y_continuous( limits = ylims, labels = scale_y) +
        scale_color_viridis_d() +
        labs(y = "Gaussian Process realisations") +
        theme_bw()

}


# SHINY START

# define server function

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
                    max = 5,
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
shinyApp(ui = ui, server = server)

