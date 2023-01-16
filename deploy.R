library(rsconnect)
deployApp(
  appFiles=c("app.R","footer.html"),
  appName = "calex",
  appTitle = "Calcium Explorer",
  logLevel = "normal"
)
