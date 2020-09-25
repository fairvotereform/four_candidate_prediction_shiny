
library(shiny)
library(tidybayes)
library(MCMCpack)
library(tidyverse)
library(cowplot)
library(glue)
library(knitr)

knit("model_explainer.Rmd", quiet = TRUE)
knit("readme.Rmd", quiet = TRUE)

ui <- fluidPage(

    # Application title
    titlePanel("Final Round RCV Poll Prediction"),

    sidebarLayout(
        sidebarPanel(
            fluidRow(
                column(5, actionButton("sim_update_button", "Simulate Poll")),
                column(5, downloadButton("download_button", 
                                      "Download Simulation"))
            ),
            tags$hr(style="border-color: black;"),
            h6("[These settings affects both prior and posterior predictions]"),
            numericInput("nA",
                         "First round count A:",
                         min = 0,
                         value = 5000),
            numericInput("nB",
                         "First round count B:",
                         min = 0,
                         value = 5000),
            numericInput("nC",
                        "First round count C:",
                        min = 0,
                        value = 3000),
            numericInput("nD",
                        "First round count D:",
                        min = 0,
                        value = 150),
            tags$hr(style="border-color: black;"),
            h6("[Settings below affect only posterior predictions]"),
            numericInput("poll_nC",
                         "Number of poll responders with first choice C:",
                         min = 0,
                         value = 50),
            numericInput("poll_nD",
                         "Number of poll responders with first choice D:",
                         min = 0,
                         value = 100),
            tags$hr(style="border-color: black;"),
            numericInput("pC_to_A",
                         "True probability that ballot from C transfers to A:",
                         min = 0,
                         max = 1,
                         value = 1/3),
            numericInput("pC_to_B",
                         "True probability that ballot from C transfers to B:",
                         min = 0,
                         max = 1,
                         value = 1/3),
            numericInput("pC_to_Exhaust",
                         "True probability that ballot from C becomes inactive:",
                         min = 0,
                         max = 1,
                         value = 1/3),
            tags$hr(style="border-color: black;"),
            numericInput("pD_to_A",
                         "True probability that ballot from D transfers to A:",
                         min = 0,
                         max = 1,
                         value = 0.1),
            numericInput("pD_to_B",
                         "True probability that ballot from D transfers to B:",
                         min = 0,
                         max = 1,
                         value = 0.2),
            numericInput("pD_to_Exhaust",
                         "True probability that ballot from D becomes inactive:",
                         min = 0,
                         max = 1,
                         value = 0.7),
            tags$hr(style="border-color: black;"),
            checkboxInput("ignoreCustomSeed", 
                          "Ignore seed input", value = TRUE),
            numericInput("seed",
                         "Seed (for random number generation): ",
                         min = 0,
                         value = as.integer(runif(1, max = 1e8))),
            textOutput("recentSeed")
            ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Posterior Predictions", 
                                 br(),
                                 h6("number of prior/posterior draws: 3,000"),
                                 hr(),
                                 h4("Simulated sample poll data"),
                                 plotOutput("sample_data"),
                                 hr(),
                                 h4("Posterior: transfer probabilities (joint)"),
                                 plotOutput("post_transfer_prob_plot"),
                                 hr(),
                                 h4("Posterior: transfer probabilities (marginal)"),
                                 plotOutput("post_transfer_prob_hist"),
                                 hr(),
                                 h4("Posterior Predictive: final round counts (joint)"),
                                 plotOutput("post_final_round_count_plot"),
                                 hr(),
                                 h4("Posterior Predictive: final round counts (marginal)"),
                                 plotOutput("post_final_round_count_hist"),
                                 hr(),
                                 h4("Posterior Predictive: final round percent"),
                                 plotOutput("post_final_round_perc_hist")),
                        tabPanel("Prior Predictions", 
                                 br(),
                                 h6("number of prior/posterior draws: 3,000"),
                                 hr(),
                                 h4("Prior: transfer probabilities (joint)"),
                                 plotOutput("prior_transfer_prob_plot"),
                                 hr(),
                                 h4("Prior: transfer probabilities (marginal)"),
                                 plotOutput("prior_transfer_prob_hist"),
                                 hr(),
                                 h4("Prior Predictive: final round counts (joint)"),
                                 plotOutput("prior_final_round_count_plot"),
                                 hr(),
                                 h4("Prior Predictive: final round counts (marginal)"),
                                 plotOutput("prior_final_round_count_hist"),
                                 hr(),
                                 h4("Prior Predictive: final round percent"),
                                 plotOutput("prior_final_round_perc_hist")),
                        tabPanel("About the model",
                                 withMathJax(includeMarkdown("model_explainer.md")))
                        )
        )
    )
)

sample_multinomial <- function(n, p1, p2, p3){
    s <- rmultinom(1, n, c(p1, p2, p3))
    data.frame(nA = s[1], nB = s[2], nExhaust = s[3])
}

sample_posterior <- function(nDraws = 3e3, 
                             nA = 100, nB = 100, nC = 100, nD = 100,
                             dirichletParamsC = c(1,1,1), 
                             dirichletParamsD = c(1,1,1)){
    cbind(
        rdirichlet(nDraws, dirichletParamsC),
        rdirichlet(nDraws, dirichletParamsD)
        ) %>%
        data.frame() %>%
        rename(pC_to_A = X1, pC_to_B = X2, pC_to_Exhaust = X3,
               pD_to_A = X4, pD_to_B = X5, pD_to_Exhaust = X6) %>%
        mutate(first_round_nA = nA,
               first_round_nB = nB,
               first_round_nC = nC,
               first_round_nD = nD,
               redistrib_C = purrr::pmap(.l = list(n = first_round_nC,
                                                    p1 = pC_to_A,
                                                    p2 = pC_to_B,
                                                    p3 = pC_to_Exhaust),
                                             .f = sample_multinomial),
               redistrib_D = purrr::pmap(.l = list(n = first_round_nD,
                                                    p1 = pD_to_A,
                                                    p2 = pD_to_B,
                                                    p3 = pD_to_Exhaust),
                                             .f = sample_multinomial)) %>%
        unnest(redistrib_C) %>%
        rename(nC_to_A = nA, nC_to_B = nB, nC_to_Exhaust = nExhaust) %>%
        unnest(redistrib_D) %>%
        rename(nD_to_A = nA, nD_to_B = nB, nD_to_Exhaust = nExhaust) %>%
        mutate(final_round_countA = first_round_nA + nC_to_A + nD_to_A,
               final_round_countB = first_round_nB + nC_to_B + nD_to_B,
               final_round_countExhaust = nC_to_Exhaust + nD_to_Exhaust,
               final_round_totalActive = final_round_countA + final_round_countB,
               final_round_percA = final_round_countA/final_round_totalActive,
               final_round_percB = final_round_countB/final_round_totalActive,
               final_round_AminusB = final_round_countA - final_round_countB,
               final_round_winnerA = final_round_countA > final_round_countB,
               draw = row_number()) %>%
        relocate(draw, .before = "pC_to_A")
        
}

plot_sample_data <- function(df, actual_params){
    
    simC <- 
        df %>%
        filter(from == "C") %>% 
        pull(count) %>%
        sum()
    
    simD <- 
        df %>%
        filter(from == "D") %>% 
        pull(count) %>%
        sum()
    
    plotC <- 
        df %>%
        filter(from == "C") %>%
        ggplot() + 
        geom_bar(aes(x = to, y = count), stat="identity") + 
        geom_point(data = actual_params %>% filter(from == "C"), 
                   aes(x = to, y = value * simC), color = "blue", size = 3) + 
        xlab("transfer from C to: ")  +
        labs(subtitle = glue("simulated final round transfer from C\n[n = {simC}]\n(blue dot is expected value)")) + 
        theme_light() + 
        theme(text = element_text(size=15))
    
    plotD <- 
        df %>%
        filter(from == "D") %>%
        ggplot() + 
        geom_bar(aes(x = to, y = count), stat="identity") + 
        geom_point(data = actual_params %>% filter(from == "D"), 
                   aes(x = to, y = value * simD), color = "blue", size = 3) + 
        xlab("transfer from D to: ")  +
        labs(subtitle = glue("simulated final round transfer from D\n[n = {simD}]\n(blue dot is expected value)")) + 
        theme_light() + 
        theme(text = element_text(size=15))
    
    plot_grid(plotC, plotD, nrow = 1)
}

plot_transfer_prob <- function(df, 
                              actual_pC_to_A, actual_pC_to_B,
                              actual_pD_to_A, actual_pD_to_B,
                              post = FALSE){
    transferC <-
        df %>%
        ggplot() + 
        geom_point(aes(x = pC_to_A, y = pC_to_B), alpha = 0.3) + 
        geom_abline(intercept = 1, slope = -1, linetype = "dashed") + 
        expand_limits(x = c(0,1), y = c(0,1)) + 
        theme_light() +
        labs(subtitle = "C transfer probability") + 
        theme(text = element_text(size=15))
    
    transferD <-
        df %>%
        ggplot() + 
        geom_point(aes(x = pD_to_A, y = pD_to_B), alpha = 0.3) + 
        geom_abline(intercept = 1, slope = -1, linetype = "dashed") + 
        expand_limits(x = c(0,1), y = c(0,1)) + 
        theme_light() +
        labs(subtitle = "D transfer probability") + 
        theme(text = element_text(size=15))
    
    if (post){
        transferC <- 
            transferC +         
            geom_point(aes(x = actual_pC_to_A, y = actual_pC_to_B), 
                       color = "blue", size = 3) +
            labs(subtitle = "C transfer probability\n(blue dot is true value)") 
            
        transferD <- 
            transferD +         
            geom_point(aes(x = actual_pD_to_A, y = actual_pD_to_B), 
                       color = "blue", size = 3) +
            labs(subtitle = "D transfer probability\n(blue dot is true value)") 
    }
    
    plot_grid(transferC, transferD, nrow = 1)
}

plot_transfer_prob_hist <- function(df, actual_params, post = FALSE){
    
    p <-
        df %>%
        pivot_longer(c(pC_to_A, pC_to_B, pC_to_Exhaust, 
                       pD_to_A, pD_to_B, pD_to_Exhaust), names_to = "param") %>%
        ggplot(aes(x = value, y = param)) + 
        stat_halfeye(point_interval = mean_qi, .width = c(0.50, 0.95)) + 
        expand_limits(x = c(0,1)) + 
        labs(subtitle = "transfer probability estimates (mean, 50%, 95% interval)") + 
        theme_light() + 
        theme(text = element_text(size=15))
    
    if (post){
        p <- 
            p + 
            geom_point(data = actual_params, aes(x = value, y = param), 
                       color = "blue", size = 3) + 
            labs(subtitle = "transfer probability estimates (mean, 50%, 95% interval)\n(blue dot is true value)") 
    }
    
    return(p)
}

plot_final_round_count <- function(df, nA, nB, nC, nD){
    df %>%
        ggplot() + 
        geom_point(aes(x = final_round_countA, y = final_round_countB), alpha = 0.3) + 
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
        expand_limits(x = c(nA, nA + nC + nD), y = c(nB, nB + nC + nD)) + 
        theme_light() +
        labs(subtitle = "predicted final round counts") + 
        theme(text = element_text(size=15))
}

plot_final_round_count_hist <- function(df, nA, nB, nC, nD){
    A <- 
        df %>%
        ggplot() + 
        geom_histogram(aes(x = final_round_countA), bins = 100,
                       color = "white", fill = "black") + 
        expand_limits(x = c(min(nA, nB), max(nA + nC + nD, nB + nC + nD))) + 
        theme_light() + 
        theme(text = element_text(size=15))
    
    B <- 
        df %>%
        ggplot() + 
        geom_histogram(aes(x = final_round_countB), bins = 100,
                       color = "white", fill = "black") + 
        expand_limits(x = c(min(nA, nB), max(nA + nC + nD, nB + nC + nD))) + 
        theme_light() + 
        theme(text = element_text(size=15))
    
    plot_grid(A, B, ncol = 1)
}


plot_final_round_perc_hist <- function(df){
    A <- 
        df %>%
        ggplot() + 
        geom_histogram(aes(x = final_round_percA), binwidth = 0.01,
                       color = "white", fill = "black") + 
        expand_limits(x = c(0,1)) + 
        labs(subtitle = glue("prob A > 50% = {mean(df$final_round_winnerA)}")) + 
        theme_light() + 
        theme(text = element_text(size=15))
    
    B <- 
        df %>%
        ggplot() + 
        geom_histogram(aes(x = final_round_percB), binwidth = 0.01,
                       color = "white", fill = "black") + 
        expand_limits(x = c(0,1)) + 
        theme_light() + 
        theme(text = element_text(size=15))
    
    plot_grid(A, B, ncol = 1)
}


# Define server logic required to draw a histogram
server <- function(input, output) {
    
    rv <- reactiveValues(seed = "none", params = NA, sim_data = NA, 
                         prior_predict = NA, post_predict = NA)

    output$sample_data <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_sample_data(rv$sim_data, rv$params)
        }
    })
        
    output$post_transfer_prob_plot <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_transfer_prob(rv$post_predict, 
                               isolate(input$pC_to_A), isolate(input$pC_to_B),
                               isolate(input$pD_to_A), isolate(input$pD_to_B),
                               post = TRUE)
        }
    })
    
    output$post_transfer_prob_hist <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_transfer_prob_hist(rv$post_predict, rv$params, post = TRUE)
        }
    })
    
    output$post_final_round_count_plot <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_final_round_count(rv$post_predict, 
                                   isolate(input$nA), isolate(input$nB),
                                   isolate(input$nC), isolate(input$nD))
        }
    })
    
    output$post_final_round_count_hist <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_final_round_count_hist(rv$post_predict, 
                                        isolate(input$nA), isolate(input$nB),
                                        isolate(input$nC), isolate(input$nD))
        }
    })
    
    output$post_final_round_perc_hist <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
           plot_final_round_perc_hist(rv$post_predict)
        }
    })
    
    output$prior_transfer_prob_plot <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_transfer_prob(rv$prior_predict, 
                               isolate(input$pC_to_A), isolate(input$pC_to_B),
                               isolate(input$pD_to_A), isolate(input$pD_to_B),
                               post = FALSE)
        }
    })
    
    output$prior_transfer_prob_hist <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_transfer_prob_hist(rv$prior_predict, rv$params, post = FALSE)
        }
    })
    
    output$prior_final_round_count_plot <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_final_round_count(rv$prior_predict, 
                                   isolate(input$nA), isolate(input$nB),
                                   isolate(input$nC), isolate(input$nD))
        }
    })
    
    output$prior_final_round_count_hist <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_final_round_count_hist(rv$prior_predict, 
                                        isolate(input$nA), isolate(input$nB),
                                        isolate(input$nC), isolate(input$nD))
        }
    })
    
    output$prior_final_round_perc_hist <- renderPlot({
        if (length(rv$sim_data) > 0 && !is.na(rv$sim_data)){
            plot_final_round_perc_hist(rv$prior_predict)
        }
    })
    
    output$recentSeed <- renderText({
        glue("most recently used seed: {rv$seed}")
    })

    observeEvent(input$sim_update_button, {
        
        if (input$ignoreCustomSeed){
            seed <- as.integer(runif(1, max = 1e8))
        }else{
            seed <- input$seed
        }
        
        rv$seed <- seed
        set.seed(seed)
        
        rv$input_df <- data.frame(
            param = c("nA", "nB", "nC", "nD",
                      "poll_nC", "poll_nD",
                      "pC_to_A", "pC_to_B", "pC_to_Exhaust",
                      "pD_to_A", "pD_to_B", "pD_to_Exhaust",
                      "seed"),
            value = c(input$nA, input$nB, input$nC, input$nD,
                      input$poll_nC, input$poll_nD, 
                      input$pC_to_A, input$pC_to_B, input$pC_to_Exhaust,
                      input$pD_to_A, input$pD_to_B, input$pD_to_Exhaust,
                      seed)
        )

        rv$params <- data.frame(param = c("pC_to_A", "pC_to_B", 
                                         "pC_to_Exhaust", "pD_to_A", 
                                         "pD_to_B", "pD_to_Exhaust"),
                                from = rep(c("C", "D"), each = 3),
                                to = rep(c("A", "B", "Exhaust"), time = 2),
                                value = c(input$pC_to_A, input$pC_to_B, 
                                         input$pC_to_Exhaust, input$pD_to_A, 
                                         input$pD_to_B, input$pD_to_Exhaust))
        
        # simulate new data
        sim_dataC <- 
            rmultinom(1, input$poll_nC, 
                      c(input$pC_to_A, input$pC_to_B, input$pC_to_Exhaust))
        
        sim_dataD <- 
            rmultinom(1, input$poll_nD, 
                      c(input$pD_to_A, input$pD_to_B, input$pD_to_Exhaust))
        
        rv$sim_data <- data.frame(from = rep(c("C", "D"), each = 3),
                                  to = rep(c("A", "B", "Exhaust"), time = 2),
                                  count = c(sim_dataC, sim_dataD))
        
        rv$prior_predict <- sample_posterior(nA = input$nA, nB = input$nB, 
                                             nC = input$nC, nD = input$nD,
                                             dirichletParamsC = c(1,1,1),
                                             dirichletParamsD = c(1,1,1))
        
        rv$post_predict <- sample_posterior(nA = input$nA, nB = input$nB, 
                                            nC = input$nC, nD = input$nD,
                                            dirichletParamsC = sim_dataC + 1,
                                            dirichletParamsD = sim_dataD + 1)
    })
    
    output$download_button <- downloadHandler(
        filename = "simulation.zip",
        content = function(file) {
            t <- tempdir()
            file.copy("readme.html", glue("{t}/readme.html"))
            setwd(t)
            write.csv(rv$input_df, "settings.csv", row.names = FALSE)
            write.csv(rv$sim_data, "sim_data.csv", row.names = FALSE)
            write.csv(rv$prior_predict, "prior.csv", row.names = FALSE)
            write.csv(rv$post_predict, "posterior.csv", row.names = FALSE)
            zip(zipfile=file, 
                files=c("readme.html", "settings.csv", "sim_data.csv", 
                        "prior.csv", "posterior.csv"))
        },
        contentType = "application/zip"
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)
