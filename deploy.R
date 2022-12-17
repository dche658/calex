library(rsconnect)
deployApp(
  appFiles=c("app.R"),
  appName = "calex",
  appTitle = "Calcium Explorer",
  logLevel = "normal"
)
