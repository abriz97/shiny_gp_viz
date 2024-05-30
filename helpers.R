library(ggplot2)
library(data.table)
library(invgamma)
library(viridis)

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
