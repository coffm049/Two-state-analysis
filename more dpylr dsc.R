library(shiny)
library(ggplot2)
library(dplyr)
library(shinydashboard)
library(nlstools)
library(zoo)
library(nls2)
df.dsc <- read.csv(file.choose(), header = TRUE)
e     <- 2.718281828459
R     <- 0.00198588


df.dsc[, seq(from = 1, to = length(df.dsc) -1, by = 2)] <- lapply(seq(from = 1, to = length(df.dsc) -1, by = 2), function(i) {
  return(if(df.dsc[1,i]>270) {df.dsc[[i]]} else {df.dsc[[i]]+273.15})
})

df.dsc[is.na(df.dsc)] <- 0


ui <- fluidPage(
  verbatimTextOutput("abla"),
  selectInput(inputId = "experiment", label ="Select dataset", choices = colnames(df.dsc)[seq(from=1, to = length(df.dsc), by =2)]),
  sliderInput(inputId = "lowtrim", label= "Select lower temp for trimming",
              min=273, max=373, value=273),
  sliderInput(inputId = "uptrim", label= "Select upper temp for trimming",
              min=273, max=373, value=373),
  downloadButton("TrimGraph Download"),
  plotOutput("trimGraph"),
  hr(),
  sliderInput(inputId = "upfbase", label= "Select upper temp for folded baseline",
              min=273, max=373, value=273),
  sliderInput(inputId = "lowubase", label= "Select lower temp for unfolded baseline",
              min=273, max=373, value=373),
  plotOutput("basegraph"),
  plotOutput("flatgraph"),
  sliderInput(inputId = "transbegin", label= "Select beginning of transition",
              min=273, max=373, value=273),
  sliderInput(inputId = "transend", label= "Select end of transition",
              min=273, max=373, value=373),
  plotOutput("geom"),
  plotOutput("nls"),
  verbatimTextOutput("conf"), 
  tableOutput("twostatetable"),
  hr(),
  checkboxInput(inputId = "hide", label= "Hide vertical lines?", value=FALSE)
)



server <- function(input, output) {
  selected  <- reactive({
    df.dsc %>%
      select(Temp = input$experiment, Cp = (which(colnames(df.dsc) == input$experiment)+1))
  })
  datasetInput <- reactive({
      switch(input$dataset,
             "rock" = rock,
             "pressure" = pressure,
             "cars" = cars)
    })
  output$trimGraph <- renderPlot({
    ggplot(data = selected(), aes(x= Temp, y=Cp)) +
      geom_point() +
      geom_vline(xintercept=input$lowtrim) +
      geom_vline(xintercept=input$uptrim) + 
      labs(x = "Temp (K)", y = "Cp(kcal/mol/K)") +
      xlim(273, 373)
  })
  trimmed <- reactive({
    selected() %>% 
      filter(Temp > input$lowtrim & Temp < input$uptrim)
  })
  fbl.Ubound.re      <- reactive({which(abs(trimmed()$Temp-input$upfbase)==min(abs(trimmed()$Temp-input$upfbase)))})
  ubl.Lbound.re      <- reactive({which(abs(trimmed()$Temp-input$lowubase)==min(abs(trimmed()$Temp-input$lowubase)))})
  output$basegraph  <- renderPlot({
    fbl             <-  lm(trimmed()$Cp[1:fbl.Ubound.re()] ~ trimmed()$Temp[1:fbl.Ubound.re()])
    ubl             <-  lm(trimmed()$Cp[ubl.Lbound.re():(length(trimmed()$Temp))] ~ trimmed()$Temp[ubl.Lbound.re():(length(trimmed()$Temp))])
    fbl.points      <- coef(fbl)[1] + coef(fbl)[2] * trimmed()$Temp
    ubl.points      <- coef(ubl)[1] + coef(ubl)[2] * trimmed()$Temp
    base.data       <- data.frame(Temp= trimmed()$Temp, fold = fbl.points, unfold = ubl.points)
    ggplot(data = trimmed(), aes(x= trimmed()$Temp, y=trimmed()$Cp)) +
      geom_point() + labs(x = "Temp (K)", y = "Cp(kcal/mol/K)") +
      geom_vline(xintercept=input$upfbase, lty=1, col="blue") +
      geom_vline(xintercept=input$lowubase, lty=1, col="blue") +
      geom_line(data=base.data, aes(x=Temp, y = fold), col="red", lty=2) +
      geom_line(data=base.data, aes(x=Temp, y = unfold), col="red", lty=2)
    })
  baselined <- reactive({
    trimmed() %>%
      mutate(Cp2 = Cp - (coef(lm(trimmed()$Cp[1:fbl.Ubound.re()] ~ trimmed()$Temp[1:fbl.Ubound.re()]))[1] + 
               coef(lm(trimmed()$Cp[1:fbl.Ubound.re()] ~ trimmed()$Temp[1:fbl.Ubound.re()]))[2] * Temp))
  })
  output$flatgraph <- renderPlot({
    ggplot(data=baselined(), aes(x=Temp, y = Cp2)) +
      geom_point()
  })
  area <- reactive({
    baselined() %>%
      filter(Temp < input$transend & Temp > input$transbegin)
  })
  dCp.geom         <- reactive({mean(tail(baselined()$Cp2,40))-mean(head(baselined()$Cp2,40))})                                 # dCp from taking difference between baselines
  int.data         <- reactive({sum(diff(area()$Temp[area()$Temp]) * rollmean(area()$Cp2[area()$Temp], 2))})
  norm.int.data    <- reactive({int.data()/tail(int.data(),n=1)})                                             # normalizes the integral to the dH of transition
  shaded           <- reactive({
    area() %>% 
      mutate(base = norm.int.data()*dCp.geom())
  })        
  trans.sig        <- reactive({area()$Cp2-shaded()$base})                    # subtracts off the spline baseline
  int.min.sig      <- reactive({cumsum(trans.sig()*diff(area()$Temp))})                   # Integrates transition without spline baseline
  norm.int.min.sig <- reactive({int.min.sig()/tail(int.min.sig(),1)})                                       # normalized integral (i.e concentration of unfolded species)
  dH.geom          <- reactive({tail(int.min.sig(), n=1) })                                               # determines dH of transion from total integral
  dH.half.geom     <- reactive({dH.geom()/2 })                                                            # determines half dH for Tm guess
  Tm.geom.index    <- reactive({which(abs(dH.half.geom()-int.min.sig())==min(abs(dH.half.geom()-int.min.sig())))})
  Tm.geom          <- reactive({area()$Temp[Tm.geom.index()]})
  H.geom           <- reactive({dH.geom()+(dCp.geom()*(baselined()$Temp-Tm.geom()))})    # define enthalpy using geom fit values
  G.geom           <- reactive({dH.geom()*(1-(baselined()$Temp/Tm.geom()))+dCp.geom()*(baselined()$Temp-Tm.geom()-baselined()$Temp*log(baselined()$Temp/Tm.geom()))})    # define free energy term
  K.geom           <- reactive({e^(-G.geom()/(R*baselined()$Temp))})                                                                  # defines  equilibrium
  dG37.geom.i      <- reactive({which(abs(baselined()$Temp-(273.15+37))==min(abs(baselined()$Temp-(273.15+37))))})
  dG37.geom        <- reactive({G.geom[dG37.geom.i()]})
  output$geom      <- renderPlot({
#    dH.df  <- data.frame(Temp=Tm.geom-7, Cp2=10)
#    dCp.df <- data.frame(Temp=Tm.geom+6, Cp2=dCp.geom/2)
#    Tm.df  <- data.frame(Temp=Tm.geom-5, Cp2=17)
#    ggplot(data=baselined(), aes(x= Temp, y = Cp2)) +
#      geom_point() + labs(x=expression(paste("Temperature ( ",degree,"C)")), y = expression(paste("Cp (kcal/mol ",degree,"C)"))) +
#      geom_point(data=shaded, aes(x=Temp, y = base), col="black") +
#      geom_ribbon(data=subset(shaded, base <= Cp2), 
#                  aes(ymin=base,ymax=Cp2), fill=rgb(0.5,0,0.1,0.25), alpha="0.5") +
#      geom_text(data=dH.df, label="Excess Enthalpy") 
    plot(baselined()$Temp-273.15, baselined()$Cp2, pch=".",xlim= c(input$lowtrim-273.15, input$uptrim-271), ylim=c(-0.1, max(baselined()$Cp2)+1), 
         xlab=expression(paste("Temperature ( ",degree,"C)")) , 
         ylab= expression(paste("Cp (kcal/mol ",degree,"C)")), cex.axis=2, cex.lab=2)
    lines(shaded()$Temp-273.15, shaded()$base)
    polygon(c(shaded()$Temp-273.15, rev(shaded()$Temp)-273.15), c(shaded()$Cp2, rev(shaded()$base)), col=rgb(0.5,0,0.1,0.25), border=NA)
    text(x=Tm.geom()-280, y=10, labels=expression(paste(Delta, "H")), cex=2, col='black')
    segments(x0=Tm.geom()-280, y0=9.5, x1=Tm.geom()-276.15, y1=4)
    segments(x0=Tm.geom()-263, y0=0, x1=Tm.geom()-263, y1=dCp.geom(), lty=3)
    text(x=Tm.geom()-267, y=dCp.geom()/2, labels=expression(paste(Delta, "Cp")), cex=2, col='black')
    segments(x0=Tm.geom()-273.15, y0=dCp.geom()/2, x1=Tm.geom()-273.15, y1=shaded()$Cp2[Tm.geom.index()], lty=2)
    text(x=Tm.geom()-5-273.15, y=17, labels="Tm", cex=2, col='black')
    segments(x0=Tm.geom()-5-273.15, y0=16.5, x1=Tm.geom()-273.15, y1=10)
    segments(x0=45, y0=0, x1=67, y1=0, lty=3, col='black')
    geom.line.fit <-  (dH.geom()*(((((dH.geom()/shaded()$Temp)-((dCp.geom()*Tm.geom())/shaded()$Temp)+dCp.geom())
                                  /(R*shaded()$Temp))*K.geom())/((K.geom()+1)^2))+((dCp.geom()*K.geom())/(K.geom()+1)))
    lines(trimmed()$Temp-273.15, geom.line.fit, col= 'red')
  })
  output$nls <- renderPlot({
    dH    <- dH.geom()
    dCp   <- dCp.geom()
    Tm    <- Tm.geom()
    fit.params <- nls(Cp2~
                        (((dH+(dCp*(Temp-Tm)))*((dH/Temp)
                                                -((dCp*Tm)/Temp)+dCp)
                          *(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                                /(0.00198588*Temp))/(0.00198588*Temp))
                          /((exp(-((dH*(1-(Temp/Tm)))
                                   +(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))/(0.00198588*Temp))+1)^2))
                         +((dCp*(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                                     /(0.00198588*Temp))))
                           /(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                                 /(0.00198588*Temp))+1))), data=baselined(),
                      start=list(Tm=Tm, dCp=dCp, dH=dH))
    ggplot(data= baselined(), aes(x= Temp, y= Cp2)) +
      geom_point() +
      geom_point(aes(x=baselined()$Temp, y = predict(fit.params)), col = 'orange', pch=".") +
      labs(x="Cp(T) (kcal/mol/K)", y="Temp(K)")
  })
  two.state.fit <- reactive({
    dH    <- dH.geom()
    dCp   <- dCp.geom()
    Tm    <- Tm.geom()
    nls(Cp2~
          (((dH+(dCp*(Temp-Tm)))*((dH/Temp)
                                  -((dCp*Tm)/Temp)+dCp)
            *(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                  /(0.00198588*Temp))/(0.00198588*Temp))
            /((exp(-((dH*(1-(Temp/Tm)))
                     +(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))/(0.00198588*Temp))+1)^2))
           +((dCp*(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                       /(0.00198588*Temp))))
             /(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                   /(0.00198588*Temp))+1))), data=baselined(),
        start=list(Tm=Tm, dCp=dCp, dH=dH))
  }) 
  output$twostatetable <- renderTable({
    dH.fit  <- summary(two.state.fit())$coefficients[3,1]
    dCp.fit <- summary(two.state.fit())$coefficients[2,1]
    Tm.fit  <- summary(two.state.fit())$coefficients[1,1]
    fitsk <- c(dH.fit, dCp.fit, Tm.fit)
    geom <- c(dH.geom(), dCp.geom(), Tm.geom())
    data.frame(Geom = geom, Fit = fitsk)
  })
  output$conf <- renderText({
    print(two.state.fit())
  })
}

shinyApp(ui = ui, server = server)

conf <-nlsConfRegions(fit.params, length=300, exp=1)
plot(conf, bounds=TRUE)
cont <- nlsContourRSS(fit.params, lseq=20, exp=5)
plot(cont, nlev=0, col=TRUE, col.pal=terrain.colors(100))