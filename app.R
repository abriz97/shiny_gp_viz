# visualise how "big" inv logit with range 1 is
library(shiny)
options(rsconnect.packrat=TRUE)
# else: remotes::install_version("rsconnect", "0.8.29")

# TODO:
# Also include plot of prior
# keep the same MVN(0, Id) draws as you vary the kernel! 
source('helpers.R')

# SHINY START

# define server function

server <- function(input, output){

    gp_realisations <- reactive({

        x <- seq(15, 50, by=1); l <- length(x)

        if(0){
            mvn_realisations <- MASS::mvrnorm(
                n = input$n_realisations,
                mu= rep(0, l),
                Sigma=diag(l) |> as.matrix()
            )
        }

        # get kernel (cholesky decompose it)
        kSE <- function(x1, x2){
            input$alpha^2 * exp( -(x1-x2)^2 / (2*input$rho^2) )
        }
        Sigma <- lapply(x, kSE, x2=x) |> 
                Reduce(f=rbind)
        if(0){
            chol_Sigma <- chol(Sigma)
            mvn_realisations <- chol_Sigma %*% t(mvn_realisations) + input$mu
        }else{
            mvn_realisations <- MASS::mvrnorm(
                n = input$n_realisations,
                mu = rep(input$mu, l),
                Sigma = Sigma 
            ) |> as.matrix()
        }
        gp_realisation <- mvn_realisations |> 
            inverse.logit() |>
            as.data.frame.matrix() |>
            as.data.table()
        names(gp_realisation) <- as.character(x)

        gp_realisation[, ID := as.factor(1:.N)]
        gp_realisation <- melt.data.table(gp_realisation,
            id.vars = "ID", variable.name="x", value.name="y") 
        gp_realisation[, x := as.numeric(x)]
    })

    plot_gp <- reactive({
        ggplot(gp_realisations(), aes(x=x, y=y, color=ID, group=ID)) + 
            geom_line() +
            geom_label(data=data.frame(),
                x=15.5, y=.95, color="black", group=NA,
                label=paste0("Mu=", input$mu)) +
            scale_y_continuous(
                limits=c(0,1), labels=scales::percent) +
            scale_color_viridis_d()+
            labs(y="Gaussian Process realisations") +
            theme_bw() 
    })

    plot_inverse_gamma <- reactive({
        x_range <- seq(0, 50, by = .05)
        data.frame(
            x=x_range,
            y = invgamma::dinvgamma(x=x_range, rate=input$rate, shape=input$shape)
        ) -> df
        # subset to not to small values
        df <- df[ df$y >= 1/10000 , ]
        ggplot(df, aes(x=x, y=y)) + 
            geom_line(color="darkred") + 
            theme_minimal() +
            labs(y="Inverse Gamma pdf")
    })

    output$gp_plot <- renderPlot({ plot_gp() })
    output$invgamma_pdf_plot <- renderPlot({ plot_inverse_gamma() })
}

# select a random bootstrap theme (for exploration)
selected_theme <- sample(bslib::bootswatch_themes(), size=1 )

# Define UI for app
ui <- shinyUI(fluidPage(

    theme = bslib::bs_theme(bootswatch=selected_theme),

    # App title ---
    titlePanel("GP hyperparameters visualisation"), 

    # Sidebar anel for inputs ---
    sidebarLayout(

        sidebarPanel( 
            width = 3, 
            conditionalPanel("input$tabs1 == 'gpreal'", 

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

            conditionalPanel("input$tabs1 == 'igprior'", 

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
            " This webpage aims to provide intuitions on how functions sampled from Gaussian Processes are affected by their hyperparameters", "\n One of the most common kernels is the squared exponential kernel, which is parametrised by the lengthscale (which I denote by the greek letter rho) and the marginal variance (denoted by alpha).\nI consider a slightly more complicated scenario whereby the function I want to model represents a probability of success, and hence has to lie in [0,1]. As such, we take an inverse logit to map the Gaussian Process samples to [0,1]. \n Happy exploration!"),

            tabsetPanel(
                type = "tabs" , id="tabs1",
                tabPanel(
                    value="gpreal",
                    id = "tab_gp_realisations",
                    title="GP realisations",
                    plotOutput(outputId="gp_plot")
                ),
                tabPanel(
                    value = "igprior",
                    id = "tab_ig_prior",
                    title="Inverse Gamma prior",
                    plotOutput(outputId="invgamma_pdf_plot")
                ),
            )
        )
    ),
    # main panel for displaying outputs

    p(em(paste("The bootstrap theme used for this webpage is:", selected_theme)))
))
shinyApp(ui=ui, server=server)

