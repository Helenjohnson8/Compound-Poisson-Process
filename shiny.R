
library(shiny)
library(ggplot2)

d_S_compound <- function(s, lambda, beta, t){
  s <- as.numeric(s)
  a <- lambda * t * beta
  out <- numeric(length(s))
  pos <- which(s>0)
  if(length(pos)>0){
    sp <- s[pos]
    out[pos] <- exp(-lambda*t - beta*sp) * sqrt(a / sp) * besselI(2*sqrt(a*sp), 1)
  }
  out[s==0] <- 0
  return(out)
}

simulate_S <- function(nsim, lambda, beta, t){
  N <- rpois(nsim, lambda*t)
  S <- numeric(nsim)
  idx.pos <- which(N>0)
  if(length(idx.pos) > 0){
    for(i in idx.pos){
      S[i] <- sum(rexp(N[i], rate=beta))
    }
  }
  S
}

ui <- fluidPage(
  titlePanel("Compound Poisson with Exponential Jumps: S(t)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("lambda", "Arrival rate λ", min=0.01, max=5, value=0.1, step=0.01),
      sliderInput("beta", "Jump rate β", min=0.1, max=10, value=1, step=0.1),
      numericInput("nsim", "Simulations (per plot)", value=20000, min=1000, max=200000, step=1000),
      checkboxInput("show_theory", "Overlay theoretical density", value=TRUE),
      selectInput("times", "Times to show", choices = c("10","100","1000","10000"),
                  selected = c("10","100","1000","10000"), multiple = TRUE),
      actionButton("go", "Simulate / Update")
    ),
    mainPanel(
      h4("Histograms (with theoretical density overlay)"),
      uiOutput("plots_ui"),
      hr(),
      h4("Summary statistics"),
      tableOutput("summary_table")
    )
  )
)

server <- function(input, output, session){
  simData <- eventReactive(input$go, {
    lambda <- input$lambda; beta <- input$beta; nsim <- input$nsim
    times <- as.numeric(input$times)
    res <- list()
    for(tt in times){
      Ssim <- simulate_S(nsim, lambda, beta, tt)
      res[[as.character(tt)]] <- Ssim
    }
    list(res=res, lambda=lambda, beta=beta, times=times)
  }, ignoreNULL = FALSE)
  
  output$plots_ui <- renderUI({
    req(simData())
    times <- simData()$times
    plot_output_list <- lapply(times, function(tt) {
      plotname <- paste0("plot_", tt)
      plotOutput(plotname, height="350px")
    })
    do.call(tagList, plot_output_list)
  })
  
  observe({
    req(simData())
    lambda <- simData()$lambda; beta <- simData()$beta
    times <- simData()$times; resList <- simData()$res
    for(tt in times){
      local({
        ttt <- tt
        plotname <- paste0("plot_", ttt)
        output[[plotname]] <- renderPlot({
          Ssim <- resList[[as.character(ttt)]]
          df <- data.frame(S=Ssim)
          maxx <- quantile(df$S, 0.995)
          p <- ggplot(df, aes(x=S)) +
            geom_histogram(aes(y=..density..), bins=80, boundary=0, closed="left", fill="grey80", color="black") +
            coord_cartesian(xlim=c(0, maxx)) +
            ggtitle(paste0("t=", ttt, "  (lambda=", lambda, ", beta=", beta, ")")) +
            theme_minimal()
          if(input$show_theory){
            sgrid <- seq(1e-8, maxx, length.out=800)
            densvals <- d_S_compound(sgrid, lambda=lambda, beta=beta, t=ttt)
            densdf <- data.frame(s=sgrid, d=densvals)
            p <- p + geom_line(data=densdf, aes(x=s, y=d), size=1)
            p <- p + annotate("text", x = maxx*0.95, y = max(densvals, na.rm=TRUE)*0.95,
                              label = paste0("P(S=0)=", round(exp(-lambda*ttt),5)),
                              hjust=1)
          }
          p
        })
      })
    }
  })
  
  output$summary_table <- renderTable({
    req(simData())
    lambda <- simData()$lambda; beta <- simData()$beta; times <- simData()$times; resList <- simData()$res
    out <- data.frame(
      t = times,
      P_S_eq_0 = sapply(times, function(tt) round(exp(-lambda*tt),6)),
      Sim_mean = sapply(times, function(tt) mean(resList[[as.character(tt)]])),
      Theor_mean = sapply(times, function(tt) round(lambda*tt/beta,6)),
      Sim_var = sapply(times, function(tt) var(resList[[as.character(tt)]])),
      Theor_var = sapply(times, function(tt) round(2*lambda*tt/beta^2,6))
    )
    out
  })
}

shinyApp(ui, server)

