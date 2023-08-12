# visualise how "big" inv logit with range 1 is
library(shiny)

# TODO:
# Also include plot of prior
# keep the same MVN(0, Id) draws as you vary the kernel!

.plot.range.alpha <- function(x_range = seq(-9, 9, by = .5), delta = 1) {
    df_base <- data.table(
        x = seq(-10, 10, by = 1 / 1000),
        y = sapply(seq(-10, 10, by = 1 / 1000), inverse.logit)
    )
    df <- data.table(
        x = x_range,
        y = sapply(x_range, inverse.logit),
        ymin = sapply(x_range - delta / 2, inverse.logit),
        ymax = sapply(x_range + delta / 2, inverse.logit)
    )
    ggplot(df_base, aes(x = x, y = y)) +
        geom_line() +
        geom_point(data = df) +
        geom_linerange(data = df, aes(ymin = ymin, ymax = ymax, color = y)) +
        theme_minimal() +
        labs(y = "inverse logit:\nexp(x)/(1+exp(x)) ")
}


# SHINY START

# define server function

server <- function(input, output) {
    require(ggplot2)
    require(data.table)
    require(invgamma)
    require(MASS)
    require(viridis)

    inverse.logit <- function(x) {
        exp(x) / (1 + exp(x))
    }

    plot.inverse.gamma <- function(rate = input$ig_rate,
                                   shape = input$ig_shape) {
        x_range <- seq(0, 50, by = .05)
        data.frame(
            x = x_range,
            y = invgamma::dinvgamma(x = x_range, rate = rate, shape = shape)
        ) -> df
        # subset to not to small values
        df <- df[df$y >= 1 / 10000, ]
        ggplot(df, aes(x = x, y = y)) +
            geom_line(color = "darkred") +
            theme_minimal() +
            labs(y = "Inverse Gamma pdf")
    }


    plot.gaussian.process.realisations <- function(mu = input$mu,
                                                   alpha = input$alpha,
                                                   rho = input$rho,
                                                   n = input$n_realisations) {
        x <- seq(15, 50, by = .5)
        l <- length(x)

        # squared exponential kernel:
        kSE <- function(x1, x2) {
            alpha^2 * exp(-(x1 - x2)^2 / (2 * rho^2))
        }
        Sigma <- lapply(x, kSE, x2 = x) |> Reduce(f = rbind)

        gp_realisation <- mvrnorm(n = n, mu = rep(0, l), Sigma = Sigma)
        gp_realisation <- gp_realisation + mu
        gp_realisation <- gp_realisation |>
            inverse.logit() |>
            as.data.frame.matrix() |>
            as.data.table()
        names(gp_realisation) <- as.character(x)

        gp_realisation[, ID := as.factor(1:.N)]
        gp_realisation <- melt.data.table(gp_realisation,
            id.vars = "ID", variable.name = "x", value.name = "y"
        )
        gp_realisation[, x := as.numeric(x)]

        ggplot(gp_realisation, aes(x = x, y = y, color = ID, group = ID)) +
            geom_line() +
            geom_label(
                data = data.frame(),
                x = 15.5, y = .95, color = "black", group = NA,
                label = paste0("Mu=", mu)
            ) +
            scale_y_continuous(
                limits = c(0, 1), labels = scales::percent
            ) +
            scale_color_viridis_d() +
            labs(y = "Gaussian Process realisations") +
            theme_bw()
    }

    output$gp_plot <- renderPlot({
        plot.gaussian.process.realisations(
            mu = input$mu,
            alpha = input$alpha,
            rho = input$rho,
            n = input$n_realisations
        )
    })
    output$invgamma_pdf_plot <- renderPlot({
        plot.inverse.gamma(rate = input$ig_rate, shape = input$ig_shape)
    })
}

# select a random bootstrap theme (for exploration)
selected_theme <- sample(bslib::bootswatch_themes(), size = 1)

# Define UI for app
ui <- shinyUI(fluidPage(
    theme = bslib::bs_theme(bootswatch = selected_theme),

    # App title ---
    titlePanel("GP hyperparameters visualisation"),

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
            p(
                " This webpage aims to provide intuitions on how functions sampled from Gaussian Processes are affected by their hyperparameters.\n One of the most common kernels is the squared exponential kernel, which is parametrised by the lengthscale (which I denote by the greek letter rho) and the marginal variance (denoted by alpha).\nI consider a slightly more complicated scenario whereby the function I want to model represents a probability of success, and hence has to lie in [0,1]. As such, we take an inverse logit to map the Gaussian Process samples to [0,1]. \n Happy exploration!"
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

