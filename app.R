require(deSolve)
require(ggplot2)
library(shiny)

ui <- fluidPage(
  sliderInput(inputId = "mu", label = "NPI efficiency", value = 1, min = 0.5, max = 1),
  sliderInput(inputId = "beta", label ="Contact rate", value = 1.1, min = 0.00, max = 10, step = 0.1),
  sliderInput(inputId = "gamma", label ="Recovery rate", value = 0.2, min = 0.00, max = 1, step = 0.1),
  sliderInput(inputId = "kappa", label ="Incubation rate", value = 0.4, min = 0.00, max = 1, step = 0.1),
  sliderInput(inputId = "f", label ="Fatigue level", value = 0.01, min = 0.0001, max = 0.1, step = 0.0001),
  
  plotOutput("plot")
)



server <- function(input, output) {
  seirmod <- function(time, variables, parameters){
    with(as.list(c(variables, parameters)), {
    lambda <- (beta*(1-mu*(0.5*tanh(x)+0.5)))
    dS <- -lambda*S*I
    dE <- lambda*S*I - kappa*E
    dI <- kappa*E - gamma*I
    dR <- gamma*I
    dx <- (-f + beta*I)
    return(list(c(dS, dE, dI, dR, dx)))
  })
  }
  
  sir_values_1 <- reactive({
    req(input$beta, input$mu, input$gamma, input$kappa, input$f)
    ode(y = c(S = 0.95, E = 0, I = 0.05, R = 0, x=0.8),
        times = seq(0, 10000, by=1),
        func = seirmod,
        parms  = c(beta=input$beta, kappa=input$kappa, gamma=input$gamma, mu=input$mu, f=input$f))
  })
  
  output$plot <- renderPlot({
    val <- as.data.frame(sir_values_1())
    
    with(val, {
      plot(time, S, type = "l", col = "blue",
           xlab = "period (days)", ylab = "number of people", ylim = c(0,1), xlim = c(0,500))
      lines(time, E, col = "orange")
      lines(time, I, col = "red")
      lines(time, R, col = "green")
    })
    
    legend("right", c("susceptibles", "exposed", "infectious", "recovered"),
           col = c("blue","orange", "red", "green"), lty = 1, bty = "n")
  })
}

shinyApp(ui = ui, server = server)

