require(ggplot2)
require(data.table)
require(invgamma)
require(MASS)
require(viridis)

inverse.logit <- function(x){
    exp(x)/(1 + exp(x))
}

.plot.range.alpha <- function(x_range = seq(-9, 9, by=.5), delta=1){
    df_base <- data.table(
        x = seq(-10,10, by=1/1000),
        y = sapply(seq(-10, 10, by=1/1000), inverse.logit)
    )
    df <- data.table(
        x = x_range,
        y = sapply(x_range, inverse.logit),
        ymin = sapply(x_range - delta/2, inverse.logit),
        ymax = sapply(x_range + delta/2, inverse.logit)
    )
    ggplot(df_base, aes(x=x, y=y)) +
        geom_line() +
        geom_point(data=df) +
        geom_linerange(data=df, aes(ymin=ymin, ymax=ymax, color=y)) + 
        theme_minimal() + 
        labs(y="inverse logit:\nexp(x)/(1+exp(x)) ")
}
