####################################
# Vinayak Rao
# For STAT695, Bayesian Data Analysis
#
# Demos for MH and Gibbs samplers
####################################

library(ggplot2)
library(shiny)

a_par <- .3
b_par <- 3

rosenbrock <- function(x,y,a=a_par,b=b_par) {
  val <- (a-x)^2 + b*(y/1-x*x)^2
  return(-val)
}

x <- seq(-2,2.5,.05)
y <- seq(-2,5.5,.05)
df <- data.frame(x=rep(x,length(y)), y = rep(y,each=length(x)))
df$z <- rosenbrock(df$x,df$y)
rosen <- ggplot(df) + geom_tile(aes(x=x,y=y,fill=exp(z))) + xlim(-2,2.5) +
  ylim(-3,6)

ii <- 1
scl <- .1
pts <- data.frame(x=rep(0,1000), y=rep(0,1000))
mcmc_type  <- 1

runApp(list(
# ui = pageWithSidebar(
  ui = fluidPage(

    headerPanel("MCMC on the Rosenbrock density"),

    selectInput("mcmc_type", label = h3("MCMC algorithm"),
                choices = list("Metropolis-Hastings" = 1, "Gibbs" = 2),
                selected = 1),

    sidebarPanel(
      sliderInput("len_scl",
                  "Lengthscale:",
                  min = .1,
                  max = 2,
                  value = .5)
    ),

    mainPanel(
      plotOutput("distPlot", click = "plot_click")
    )
  ),
  server =function(input, output, session) {
    autoInvalidate <- reactiveTimer(500, session)

    output$distPlot <- renderPlot({
      autoInvalidate()
      # generate an rnorm distribution and plot it
      ii       <<- ii + 1
      if(mcmc_type == 1) {
        prop     <- pts[ii-1,] + scl * rnorm(2);
        pts[ii,] <<- prop
        if(log(runif(1)) > rosenbrock(pts[ii,1], pts[ii,2]) - rosenbrock(pts[ii-1,1],pts[ii-1,2]))
          pts[ii,] <<- pts[ii-1,]
        rosen + geom_path(data = pts[1:ii,], aes(x=x,y=y), color='red',size=2, alpha=.4) +
          geom_point(data=prop, aes(x=x,y=y), size=4,color='black')
      } else {
        pts[ii,] <<- pts[ii-1,]
        if(ii%%2) {
          fx     <- function(x) {exp(rosenbrock(x,pts[ii,2]))}
          z      <- integrate(fx,-Inf,Inf)$value
          fx_cdf <- function(x) {integrate(fx,-Inf,x)$value/z - runif(1)}
          pts[ii,1] <<- uniroot(fx_cdf,lower = -6, upper = 6)$root
        } else {
          pts[ii,2] <<- rnorm(1, pts[ii,1]^2, .5/sqrt(b_par))
        }
        rosen + geom_path(data = pts[1:ii,], aes(x=x,y=y), color='red',size=2, alpha=.4) +
          geom_point(data=pts[ii,], aes(x=x,y=y), size=4,color='black')
      }
    })

   observeEvent(input$plot_click, {
     pts[1,] <<- c(input$plot_click$x, input$plot_click$y)
     ii      <<- 1
   })

    observeEvent(input$len_scl, {
      scl <<- input$len_scl
    })

    observeEvent(input$mcmc_type, {
      mcmc_type <<- input$mcmc_type
    })

 }
))
