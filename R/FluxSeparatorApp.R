#' @title FluxSeparator Shiny App
#'
#' @description
#' Shiny app to visualize the process of separating diffusive and ebullitive fluxes. The app also gives the user the ability to
#' mark areas where diffusive fluxes should be calculated.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{shinyApp}()}.
#'
#' @references Sø et al. (2024). Self-Made Equipment for Automatic Methane
#'   Diffusion and Ebullition Measurements From Aquatic Environments.
#'   \doi{10.1029/2024JG008035}.
#' @references Sø et al. (2023). Methane and carbon dioxide fluxes at high
#'   spatiotemporal resolution from a small temperate lake.
#'   \doi{10.1016/j.scitotenv.2023.162895}.
#' @references Hutchinson, G.L. and Mosier, A.R. (1981). Improved soil cover
#'   method for field measurement of nitrous oxide fluxes. Soil Science Society
#'   of America Journal, 45, pp. 311-316.
#' @references Pullens, J.W.M., Abalos, D., Petersen, S.O. and Pedersen, A.R.
#'   (2023). Identifying criteria for greenhouse gas flux estimation with
#'   automatic and manual chambers: A case study for N2O. European Journal of
#'   Soil Science, 74, e13340. \doi{10.1111/ejss.13340}.
#'
#' @seealso \code{\link{diffusive_flux}}, \code{\link{ebullitive_flux}}, \code{\link{ppm_to_umol}}
#'
#' @author Jonas Stage Sø \email{Jonassoe@biology.sdu.dk}
#'
#' @importFrom rlang .data
#' @import shiny
#' @import shinydashboard
#' @importFrom DT DTOutput renderDT
#' @importFrom dplyr rename_with distinct count case_when lead where everything starts_with rename
#' @importFrom stringr str_to_lower
#' @importFrom lubridate ymd_hm ydm_hms mdy_hms dmy_hms
#' @importFrom readr read_delim col_integer cols
#' @importFrom ggplot2 geom_line theme_bw sec_axis theme
#' @importFrom stats coef na.exclude var
#' @importFrom utils data head write.csv
#'
#' @examples
#' \donttest{
#' library(FluxSeparator)
#'
#' FluxSeparatorApp()
#' }
#' @export

FluxSeparatorApp <- function(...) {

  ggplot2::theme_set(theme(panel.background = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           legend.text = element_text(size = 14),
                           legend.title = element_text(size = 14),
                           legend.position = "bottom",
                           plot.title = element_text(hjust = 0.5),
                           axis.text = element_text(size = 12),
                           axis.title = element_text(size = 14),
                           strip.background = element_rect(fill = "white"),
                           panel.border = element_rect(colour = "black", fill=NA),
                           legend.key = element_rect(fill = "NA",color = "white"),
                           legend.background = element_rect(color = "white"),
                           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), angle = 90),
                           axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)),
                           strip.placement = "outside",
                           strip.text = element_text(size = 14),
                           strip.text.x = element_blank()))

  options(shiny.maxRequestSize = 100*1024^2)

  read_csv(system.file("shiny/Data/all_sensor_model_coef.csv", package = "FluxSeparator")) %>%
    rename_with(str_to_lower,everything()) -> model_coef

# UI ----
  ui <- dashboardPage(
    dashboardHeader(title = "Flux calculations"),
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        menuItem("Welcome", tabName = "intro", icon = icon("door-open")),
        menuItem("Start", tabName = "start", icon = icon("play")),
        menuItem("Ebullitive flux", tabName = "ebul", icon = icon("circle")),
        menuItem("Diffusive flux", tabName = "diff", icon = icon("chart-line"))
      )
    ),
    dashboardBody(
      tabItems(
        # Intro UI tab ----
        tabItem(tabName = "intro",
                HTML("
                     <h2>Introduction to the flux calculator</h2>
                     Welcome to this website which helps to calculate fluxes.
                     The website allows the user to separate ebullitive (bubble) and diffusive fluxes.
                     The technique for separation ebullitive and diffusive fluxes is based on the R-package <a href='https://github.com/JonasStage/FluxSeparator'> <em>FluxSeparator</em></a> with minor alterations.
                     Similar results can be obtained by using the R-package in R, however, this website allows for easy usage and testing how different parameters affect your fluxes.<br>

                     This website was originally made to ease the usage of DIY sensors for <a href='https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024JG008035'> CH<sub>4</sub> and CO<sub>2</sub></a> supplied by the author.
                     Nevertheless, increasing interest fueled me to incorporate other types of sensors.<br>

                     <h2>Where to start</h2>
                     Start by selecting the 'Start' tab on the left side, here you can upload your data and make sure the formatting is correct.
                     Before uploading the data a selection can be made whether to use 'DIY' sensors, which refers to paper just above, or 'Other' which can be any format.
                     The 'Start' page also displays your data, allowing you to zoom in on particular areas and calculate the flux from these areas.
                     After uploading the data and ensuring it all looks good, you can head on to one of the next tabs for calculating either ebullitive or diffusive fluxes. <br>

                     <h2>Theory behind the calculations </h2>
                     <h3>Ebullitive flux</h3>
                     Most notable is the runvar_cutoff variable which determines the threshold value of the running variance, with running variances above this threshold being considered bubbles.
                     Furthermore, selections allow the user to include more or less observations before and after the bubble, or set a minimum threshold for the concentration change caused by the ebullitive event.<br>
                     <h3>Diffusive flux</h3>
                     The diffusive flux allows the user to also consider the bubbles from the ebullitive flux function.
                     If this is selected the function will only calculate diffuisve fluxes before any bubbles occur, to ensure concentrations are not elevated.
                     A possibility to not look for bubbles is also present, which is useful when looking at CO<sub>2</sub> fluxes.
                     For more information on how to use the different parameters and how the separation of fluxes is done, see <a href='https://doi.org/10.1029/2024JG008035'> this paper</a> <br><br>

                     <h2>Citing this tool</h2>
                     Please cite the <a href='https://doi.org/10.1029/2024JG008035'> paper introducing the R-package</a> and the <a href='https://doi.org/10.5281/zenodo.8297154'> R-package </a> when using this tool: <br><br>
                     <em>S\u00f8, J. S., Sand-Jensen, K., & Kragh, T. (2024). Self-made equipment for automatic methane diffusion and ebullition measurements from aquatic environments. Journal of Geophysical Research: Biogeosciences, 129, e2024JG008035. https://doi.org/10.1029/2024JG008035</em><br><br>
                     <em>S\u00f8, J. S., Sand-Jensen, K., & Kragh, T. (2023). FluxSeparator - Separation of diffusive and ebullitive fluxes (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.8297154</em><br>

                     <h2>Wrapping up</h2>

                     I expect there will be encounters of errors while using the website, which I would gladly try to accomodate, so please let me know.
                     Additionally, don't hesitate to contact me if you have ideas on how to improve this website.<br><br>
                     This website is made by <a href='mailto:Jonassoe@biology.sdu.dk'</a> Jonas Stage S\u00f8 </a>, Ph.D., Postdoc at The University of Southern Denmark. <br><br>

                     <center>
                     <a href='https://www.researchgate.net/profile/Jonas-Stage-So?ev=hdr_xprf'><i class='fa-brands fa-researchgate' style='font-size: 6em'></i></a>
                     <a href='https://github.com/JonasStage'><i class='fa-brands fa-github' style='font-size: 6em'></i></a><br><br><br>
                     </center>

                     <h6><em>We thank Kenneth T. Martinsen for help developing this app</em></h6>"),
        ),
        # Start UI tab ----
        tabItem(tabName = "start",
                titlePanel("Flux sensor calculations"),

                sidebarLayout(

                  sidebarPanel(
                    radioButtons("file_type", "What type of instrument-file is being used",
                                 choices = c("DIY","Other")),

                    tags$b("Upload data file"),

                    conditionalPanel(
                      condition = "input.file_type == 'DIY'",

                      fileInput("file", "",
                                multiple = FALSE,
                                accept = c("text/csv",
                                           "text/comma-separated-values",
                                           ".csv")),

                      radioButtons("calibration_choice","Use own or Jonas' calibration file",
                                   choices = c("Jonas","Own")),

                      conditionalPanel(
                        condition = "input.calibration_choice == 'Jonas'",
                        textInput("sensor_name","Input sensor name for calibration")),

                      conditionalPanel(
                        condition = "input.calibration_choice == 'Own'",
                        tags$b("Upload own calibration file"),
                        tags$p("The file should be in csv format with a header and should include the following columns: a, b, c, k, g and s"),
                        fileInput("own_sensor_cal_file",
                                  label = NULL,
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values",
                                             ".csv"))),

                      tags$b("Timezone of measurement (UTC)"), numericInput("timezone",
                                                                            label = "",
                                                                            value = 0,
                                                                            min = -12,
                                                                            max = 14),
                    ),
                    conditionalPanel(
                      condition = "input.file_type == 'Other'",
                      fileInput("file_other", "",
                                multiple = FALSE),
                      radioButtons("file_format", "Select file format",
                                   choices = c(".csv",
                                               ".txt",
                                               ".data")),
                      tags$b("How many rows to skip when reading in data"), numericInput("skip_rows",
                                                                                         label = "",
                                                                                         value = 0,
                                                                                         min = 0,
                                                                                         max = Inf),
                      radioButtons("datetime_format","Select datetime format",
                                   choices = c("YMD-HMS",
                                               "YDM-HMS",
                                               "MDY-HMS",
                                               "DMY-HMS")),
                      selectInput("datetime_column","Select datetime column",
                                  choices = c(colnames(data()))),
                      selectInput("concentration_values_ch4_column","Select first concentration column",
                                  choices = c(colnames(data()))),
                      radioButtons("unit_concentration1", "Select units for the first concentration column",
                                   choices = c("ppm","ppb")),
                      selectInput("airt_column","Select temperature column (\u00b0C)",
                                  choices = c(colnames(data()),NA_character_)),
                      selectInput("concentration_values_co2_column","Select second concentration column if applicable",
                                  choices = c(colnames(data()),NA_character_)),
                      radioButtons("unit_concentration2", "Select units for the second concentration column",
                                   choices = c("ppm","ppb")),
                      selectInput("water_column","Select water vapor concentration column if applicable",
                                  choices = c(colnames(data()),NA_character_)),
                      selectInput("sep_column","Select measurement separator column if applicable",
                                  choices = c(colnames(data()),NA_character_))
                    ),
                    tags$b("Download data as csv"),

                    tags$br(),

                    #Button to save the file
                    downloadButton('ch4_download', 'Download'),

                    tags$hr(style = "border-top: 3px solid #000000;"),

                    tags$b("Plotting options"),

                    tags$p("Adjust time range"),

                    sliderInput("range", "",
                                ymd_hm("2024-01-01 12:00"),
                                ymd_hm("2024-12-31 12:00"),
                                c(ymd_hm("2024-01-01 12:00"),
                                  ymd_hm("2024-12-31 12:00")),
                                60*60*24,
                                timezone="+0000"),

                    tags$p("Adjust first concentration range"),

                    sliderInput("ch4_range", "",
                                min = 0, max = 10000, value = c(0, 10000)),

                    tags$hr(),

                    tags$p(HTML(paste0("Adjust second contration range"))),

                    sliderInput("co2_range", "",
                                min = 0, max = 10000, value = c(0, 10000)),

                    tags$hr(style = "border-top: 3px solid #000000;"),

                    tags$b("Enter chamber metadata"),
                    fluidRow(
                      column(4,
                             tags$p("Chamber volume (L):")),
                      column(4,
                             tags$p(HTML(paste0("Chamber area (m",tags$sup(2),"):")))),
                      column(4,
                             tags$p("Atmos. pressure (atm):"))),

                    fluidRow(
                      column(4,
                             numericInput("chamber_vol", NULL, 10, min = 0)),
                      column(4,
                             numericInput("chamber_area", NULL, 0.5, min = 0)),
                      column(4,
                             numericInput("atm_pres", NULL, 1, min = 0))),

                    textInput("first_con_name", "First concentration name"),
                    textInput("second_con_name", "Second concentration name"),

                  ),

                  mainPanel(
                    conditionalPanel(
                      condition = "input.file_type == 'Other'",
                      h5("A preview of you data will show here to help ensure data is correctly formatted. If the preview does not appear as a table try altering the amount of rows that are to be skipped"),
                      tableOutput("head_df")
                    ),

                    tags$b("Main plot"),
                    p("Use sliders to adjust the x-axis"),

                    plotOutput("plot",
                               brush = brushOpts(
                                 id = "plot_brush",
                                 resetOnNew = TRUE
                               )),

                    tags$hr(style = "border-top: 3px solid #000000;"),

                    tags$b("Zoom plot"),
                    p("Use mouse to select measurements in the main plot"),

                    plotOutput("plot_zoom"),

                    htmlOutput("result_string"),

                    tags$hr(style = "border-top: 3px solid #000000;"),

                    tags$b("Save flux of zoom plot"),

                    tags$p("Flux ID (optional):"),

                    fluidRow(
                      column(4,
                             textInput("sample_id", NULL, "ID", width = "200px")),
                      column(2,
                             actionButton("save", "Save")),
                      column(2,
                             actionButton("delete_saved", "Delete all saved data"))),

                    tags$b("Saved data"),
                    p("Table with saved data, export table as '.csv' file using the download button"),
                    DTOutput("results"),

                    tags$b("Download flux data as csv"),

                    tags$br(),

                    #Button to save the file
                    downloadButton('download', 'Download'),

                    tags$hr(style = "border-top: 3px solid #000000;"),
                  ))),
        # Ebul UI tab ----
        tabItem(tabName = "ebul",
                titlePanel("Ebullitive flux"),
                sidebarLayout(
                  sidebarPanel(
                    selectInput("concentration_values","Select concentration column",
                                choices = c(colnames(data())),
                                selected = "ch4"),
                    sliderInput("runvar_cutoff","Running variance cutoff",
                                min = 0 , max = 2, value = 0.5,step = 0.001),
                    selectInput("top_selection", "Select how to calculate the ebullitive flux",
                                choices = c("Last" = "last","Maximum" = "max")),
                    sliderInput("indexspan", "Decide on number of observations before and after the running variance cutoff",
                                min = 0, max = 100, value = 30, step = 1),
                    sliderInput("concentration_diffusion_cutoff", "Decide on the cutoff value for diffusion",
                                min = 0, max = 100, value = 1, step = 1),
                    sliderInput("number_of_pumpcycles_in_plot", "How many plots to show at the same time",
                                min = 1, max = 24, value = 10, step = 1),
                    tags$b("Choose whether or not to smooth data"),
                    checkboxInput("smooth_data","Smooth data",value = F),

                  ),
                  mainPanel(
                    h3("Start by uploading your data in the Start tab"),
                    plotOutput("plot2"),
                    column(
                      12,
                      sliderInput("plot_number_select", "Switch between plots",
                                  min = 1, max = 100, value = 1, step = 1),
                      align = "center"),
                    DTOutput("ebul_table"),
                    column(
                      12,
                      downloadButton('download_ebul', 'Download shown fluxes'),
                      align = "center")
                  ))),
        # Diff UI tab ----
        tabItem(tabName = "diff",
                titlePanel("Diffusive flux"),
                sidebarLayout(
                  sidebarPanel(
                    selectInput("concentration_values_diff","Select concentration column",
                                choices = c(colnames(data())),
                                selected = "CH4"),
                    tags$b("Only calculate diffusive fluxes where there has been no previous bubbles"),
                    checkboxInput("look_for_bubbles","Look for bubbles",value = T),
                    sliderInput("remove_observations_prior", "How many observations to remove prior to calculations of diffusive flux",
                                min = 0, max = 10000, value = 1, step = 1),
                    sliderInput("number_of_observations_used", "How many observations to use for the calculations of diffusive flux",
                                min = 0, max = 10000, value = 400, step = 1),
                    sliderInput("number_of_observations_required", "How many are required before diffusive flux is calculated",
                                min = 0, max = 10000, value = 50, step = 1),
                    sliderInput("cutoff_start_value", "The upper value that a starting concentration can be",
                                min = 0, max = 10000, value = 100, step = 1),
                    sliderInput("number_of_pumpcycles_in_plot_diff", "How many plots to show at the same time",
                                min = 1, max = 24, value = 10, step = 1),
                    tags$b("Choose whether or not to smooth data"),
                    checkboxInput("smooth_data_diff","Smooth data",value = F),
                    tags$b("Choose whether or not to apply the Hutchinson-Mosier correction"),
                    checkboxInput("hmr_correction","Hutchinson-Mosier correction",value = F),
                    tags$b("Supply volume and area if Hutchinson-Mosier correction is to be applied"),
                    fluidRow(
                      column(
                        6,numericInput("volume_diff", "Volume (L)", 10, min = 0)
                      ),
                      column(
                        6, numericInput("area_diff", HTML(paste0("Chamber area (m",tags$sup(2),"):")), 0.5, min = 0),
                      ))),
                  mainPanel(
                    h3("To calculate diffusive fluxes, ensure you have also considered the ebullitive fluxes on the previous page or select to not look for bubbles"),
                    plotOutput("plot_diff"),
                    column(
                      12,
                      sliderInput("plot_number_select_diff", "Switch between plots",
                                  min = 1, max = 100, value = 1, step = 1),
                      align = "center"),
                    DTOutput("diff_table"),
                    column(
                      12,
                      downloadButton('download_diff', 'Download shown fluxes'),
                      align = "center"),
                  )))
      )))

# Server ----
  server <- function(input, output, session){
    data <- reactive({

      if(input$file_type == "DIY"){
        req(input$file)
        req(input$timezone)
        updateTextInput(session, "first_con_name", value = "CH4")
        updateTextInput(session, "second_con_name", value = "CO2")
        if(input$calibration_choice == "Jonas"){
          req(input$sensor_name)
          model_coef %>%
            mutate(sensor = str_to_lower(sensor)) %>%
            filter(sensor == str_to_lower(input$sensor_name)) -> calibration_constants
          validate(
            need(dim(calibration_constants)[1] != 0,
                 message = "Unrecognized sensor name")
          )
        } else {
          req(input$own_sensor_cal_file)
          read_csv(input$own_sensor_cal_file$datapath) %>%
            rename_with(str_to_lower,everything()) -> calibration_constants
          validate(
            need("a" %in% colnames(calibration_constants),
                 message = "The 'a' value is missing in the uploaded calibration file"),
            need("b" %in% colnames(calibration_constants),
                 message = "The 'b' value is missing in the uploaded calibration file"),
            need("c" %in% colnames(calibration_constants),
                 message = "The 'c' value is missing in the uploaded calibration file"),
            need("k" %in% colnames(calibration_constants),
                 message = "The 'k' value is missing in the uploaded calibration file"),
            need("g" %in% colnames(calibration_constants),
                 message = "The 'g' value is missing in the uploaded calibration file"),
            need("s" %in% colnames(calibration_constants),
                 message = "The 's' value is missing in the uploaded calibration file"),
          )
        }
        req(calibration_constants)
        lookup <- c(rh = "RH")

        df <- read_csv(input$file$datapath, col_names = F,
                       col_types = cols(X1 = col_integer(),
                                        X2 = col_integer(),
                                        X3 = col_character(),
                                        X4 = col_double(),
                                        X5 = col_double(),
                                        X6 = col_double(),
                                        X7 = col_double(),
                                        X8 = col_double(),
                                        X9 = col_double(),
                                        X10 = col_double(),
                                        X11 = col_double(),
                                        X12 = col_integer(),
                                        X13 = col_integer())) %>%
          rename(millis = X1,
                 stampunix = X2,
                 datetime = X3,
                 RH = X4,
                 tempC = X5,
                 CH4smV = X6,
                 CH4rmV = X7,
                 VbatmV = X8,
                 K33_RH = X9,
                 K33_Temp = X10,
                 K33_CO2 = X11,
                 SampleNumber = X12,
                 PumpCycle = X13) %>%
          rename(rh = starts_with("RH"),
                 ch4_smv=CH4smV) %>%
          cbind(calibration_constants) %>%
          filter(!is.na(datetime),
                 between(rh, 0,100)) %>%
          filter(lead(!is.na(SampleNumber)), !is.na(SampleNumber)) %>%
          mutate(datetime = ymd_hms(datetime),
                 airt = as.numeric(tempC),
                 abs_H = (6.112*exp((17.67*airt)/(airt+243.5))*rh*18.02)/((273.15+airt)*100*0.08314),
                 ppm_H20 = 1358.326542*abs_H,
                 second_concentration = (K33_CO2/(1-(ppm_H20/10^6))),
                 V0 = abs_H*g+s,
                 RsR0 = ((5000/ch4_smv)-1)/((5000/V0)-1),
                 first_concentration = a*(RsR0^b)+c*abs_H*(a*RsR0^b) + k,
                 datetime = datetime+(input$timezone-1)*3600) %>%
          rename(water = ppm_H20) %>%
          group_by(datetime) %>%
          summarise(across(c(rh, airt, second_concentration, first_concentration, water, PumpCycle), mean)) %>%
          select(datetime, rh, airt, second_concentration, first_concentration, water,PumpCycle)
      } else {
        file_upload <- raw_file_other()

        if(!"datetime" %in% colnames(file_upload) & ncol(file_upload) > 6 & "date" %in% colnames(file_upload) & "time" %in% colnames(file_upload)){
          file_upload <- file_upload %>%
            mutate(datetime = paste(date,time))
        }

        datetime_format <- input$datetime_format
        unit1 <- input$unit_concentration1
        unit2 <- input$unit_concentration2

        file_upload %>%
          cbind(datetime_format,unit1,unit2) %>%
          rename(datetime = any_of(input$datetime_column),
                 first_concentration = any_of(input$concentration_values_ch4_column),
                 second_concentration = any_of(input$concentration_values_co2_column),
                 water = any_of(input$water_column),
                 airt = any_of(input$airt_column),
                 PumpCycle = any_of(input$sep_column)) %>%
          mutate(datetime = case_when(datetime_format == "YMD-HMS" ~ ymd_hms(datetime),
                                      datetime_format == "YDM-HMS" ~ ydm_hms(datetime),
                                      datetime_format == "MDY-HMS" ~ mdy_hms(datetime),
                                      datetime_format == "DMY-HMS" ~ dmy_hms(datetime)),
                 across(any_of(c("first_concentration","second_concentration","water","airt","PumpCycle")), ~ parse_number(.x))) %>%
          select(any_of(c('datetime', 'airt', 'second_concentration',"first_concentration","water","PumpCycle"))) -> df


        if(input$unit_concentration1 == "ppb") {
          df <- df %>%
            mutate(first_concentration = first_concentration/1000)
        }
        if(input$unit_concentration2 == "ppb") {
          df <- df %>%
            mutate(second_concentration = second_concentration/1000)
        }


        if(!"PumpCycle" %in% colnames(df)) {
          df <- mutate(df, PumpCycle = 1)
        }
        if(!"co2" %in% colnames(df)){
          df <- mutate(df, co2 = 0)
        }
        if(!"water" %in% colnames(df)){
          df <- mutate(df, water = 0)
        }
        if(!"airt" %in% colnames(df)){
          df <- mutate(df, airt = 0)
        }
      }

      time_start <- min(df$datetime, na.rm =T)
      time_end <- max(df$datetime, na.rm =T)

      ch4_start <- floor(min(df$first_concentration, na.rm=T))
      ch4_end <- ceiling(max(df$first_concentration, na.rm=T))

      co2_start <- floor(min(df$second_concentration, na.rm =T))
      co2_end <- ceiling(max(df$second_concentration, na.rm =T))

      updateSliderInput(session, "range", value = c(time_start, time_end),
                        min = time_start, max = time_end, step = 60)

      updateSliderInput(session, "ch4_range", value = c(ch4_start, ch4_end),
                        min = ch4_start, max = ch4_end, step = 1)

      updateSliderInput(session, "co2_range", value = c(co2_start, co2_end),
                        min = co2_start, max = co2_end, step = 1)

      if(nchar(input$first_con_name)>0) {
        names(df)[names(df) == "first_concentration"] <- input$first_con_name
      }

      if(nchar(input$second_con_name)>0) {
        names(df)[names(df) == "second_concentration"] <- input$second_con_name
      }

      updateSelectInput(session, "concentration_values", choices = colnames(df), selected = "CH4")
      updateSelectInput(session, "concentration_values_diff", choices = colnames(df), selected = "CH4")
      updateSelectInput(session, "datetime_column", choices = colnames(df), selected = "NA_character_")

      return(df)
    })

    to_listen<- reactive({
      list(input$file_other,
           input$skip_rows,
           input$file_format)
    })

    raw_file_other <- reactive({
      req(input$file_other)
      req(input$file_format)
      if (input$file_format == ".csv") {
        read_csv(input$file_other$datapath,
                 skip = input$skip_rows,
                 col_types = cols(.default = col_character()))
      } else if (input$file_format %in% c(".txt", ".data")) {
        read_delim(input$file_other$datapath,
                   delim = "\t",
                   skip = input$skip_rows,
                   col_types = cols(.default = col_character()))
      }
    })

    observeEvent(to_listen(), {
      req(raw_file_other())
      file_upload <- raw_file_other()

      updateSelectInput(session, "datetime_column", choices = c(colnames(file_upload),"datetime"), selected = "NA_character_")
      updateSelectInput(session, "concentration_values_ch4_column", choices = c(colnames(file_upload),NA_character_), selected = "NA_character_")
      updateSelectInput(session, "concentration_values_co2_column", choices = c(colnames(file_upload),NA_character_), selected = "NA_character_")
      updateSelectInput(session, "water_column", choices = c(colnames(file_upload),NA_character_), selected = "NA_character_")
      updateSelectInput(session, "airt_column", choices = c(colnames(file_upload),NA_character_), selected = "NA_character_")
      updateSelectInput(session, "sep_column", choices = c(colnames(file_upload),NA_character_), selected = "NA_character_")

      output$head_df<- renderTable({
        file_upload %>%
          head(5)
      })
    })

    observeEvent(input$concentration_values_co2_column, {
      updateTextInput(session, "second_con_name", value = input$concentration_values_co2_column)
    }, ignoreInit = TRUE)

    data_out <- reactiveVal(data.frame())

    output$plot <- renderPlot({
      if(input$file_type == "Other") {
        validate(
          need(input$datetime_column, message = "Needs to input datetime column"),
          need(input$concentration_values_ch4_column, message = "Needs to input first concentration column"))

        updateTextInput(session, "first_con_name", value = input$concentration_values_ch4_column)
      }

      req(data())
      req(input$first_con_name)
      data_subset <- data() %>%
        filter(between(datetime, input$range[1], input$range[2]),
               between(.data[[input$first_con_name]], input$ch4_range[1],input$ch4_range[2]))

      if(nchar(input$second_con_name)>0){
        data_subset <- data_subset %>%
          filter(between(.data[[input$second_con_name]], input$co2_range[1],input$co2_range[2]))
      }

      ggplot() +
        geom_point(data = data_subset, aes(.data$datetime, .data[[input$first_con_name]], col = "CH4")) +
        labs(x = "Datetime",
             y = input$first_con_name,
             col = "") +
        scale_color_manual(
          labels = input$first_con_name,
          values = c("darkorange")) +
        scale_x_datetime(date_breaks = "10 min",
                         date_minor_breaks = "1 min",
                         date_labels = "%R") +
        theme_bw() +
        theme(legend.position = "bottom") -> op1

      if (nchar(input$second_con_name)==0) {
        op1
      } else {

        co2_min = min(data_subset[[input$second_con_name]], na.rm=T)
        co2_max = max(data_subset[[input$second_con_name]], na.rm=T)
        ch4_min = min(data_subset[[input$first_con_name]], na.rm=T)
        ch4_max = max(data_subset[[input$first_con_name]], na.rm=T)

        ch4_scaled = (co2_max - co2_min)*((data_subset[[input$first_con_name]]-ch4_min)/(ch4_max - ch4_min))+co2_min
        ch4_labels = pretty(data_subset[[input$first_con_name]])
        ch4_at = (co2_max - co2_min)*((ch4_labels-ch4_min)/(ch4_max - ch4_min))+co2_min

        co2_scaled = (ch4_max - ch4_min)*((data_subset[[input$second_con_name]]-co2_min)/(co2_max - co2_min))+ch4_min
        co2_labels = pretty(data_subset[[input$second_con_name]])
        co2_at = (ch4_max - ch4_min)*((co2_labels-co2_min)/(co2_max - co2_min))+ch4_min

        op1 +
          geom_point(data = data_subset, aes(datetime, co2_scaled, col = "CO2")) +
          labs(x = "Datetime",
               y = input$first_con_name,
               col = "") +
          scale_y_continuous(sec.axis = sec_axis(trans=~., name = input$second_con_name,
                                                 breaks = co2_at, labels = co2_labels)) +
          scale_color_manual(limits = c("CH4","CO2"),
                             labels = c(input$first_con_name, input$second_con_name),
                             values = c("darkorange", "forestgreen")) +
          scale_x_datetime(date_breaks = "1 hour",
                           date_minor_breaks = "30 min",
                           date_labels = "%R")
      }

    })

    ranges2 <- reactiveValues(x = NULL, y = NULL)

    data_subset <- reactive({
      if(input$file_type == "DIY") req(input$file) else req(input$file_other)

      if (!is.null(ranges2$x)) {
        ranges2$x <- as_datetime(ranges2$x)

        data_subset <- data() %>%
          filter(between(datetime, ranges2$x[1], ranges2$x[2])) %>%
          mutate(sec = cumsum(c(0, diff(as.numeric(datetime)))))

        lm_model_ch4 <- lm(data_subset[[input$first_con_name]]~sec, data = data_subset)
        slope_ch4 <- coef(lm_model_ch4)[2]
        intercept_ch4 <- coef(lm_model_ch4)[1]
        r2_ch4 <- summary(lm_model_ch4)$r.squared

        mean_temp <- mean(data_subset$airt)

        R <- 0.08206 #L atm K^-1 mol^-1


        ch4_flux <- (slope_ch4*(input$chamber_vol*input$atm_pres))/(R*(273.15+mean_temp)*input$chamber_area)

        if(var(data_subset$airt) == 0) {
          ch4_flux = NA_real_
        } else {}

        if (nchar(input$second_con_name)==0) {

          results_string <- paste0(input$first_con_name,": slope = ", round(slope_ch4*3600, 2), " (ppm h<sup>-1</sup>)",
                                   ", flux = ", round(ch4_flux*3600, 2), " (\u03bcmol m<sup>-2</sup> h<sup>-1</sup>)",
                                   ", R<sup>2</sup> = ", round(r2_ch4, 2))

          results <- data.frame("processing_date" = strftime(Sys.time(), "%Y-%m-%d %H:%M:%S", tz="GMT"),
                                "id" = as.character(input$sample_id),
                                "start" = strftime(ranges2$x[1], "%Y-%m-%d %H:%M:%S", tz="GMT"),
                                "end" = strftime(ranges2$x[2], "%Y-%m-%d %H:%M:%S", tz="GMT"),
                                "concentration1_name" = input$first_con_name,
                                "concentration1_slope_h1" = slope_ch4*3600,
                                "concentration1_intercept" = intercept_ch4,
                                "concentration1_R2" = r2_ch4,
                                "temperature" = mean_temp,
                                "chamber_volume" = input$chamber_vol,
                                "chamber_area" = input$chamber_area,
                                "concentration1_flux_umol_m2_h1" = ch4_flux*3600)

        } else {
          lm_model_co2 <- lm(data_subset[[input$second_con_name]]~sec, data = data_subset)
          slope_co2 <- coef(lm_model_co2)[2]
          intercept_co2 <- coef(lm_model_co2)[1]
          r2_co2 <- summary(lm_model_co2)$r.squared

          co2_flux <- (slope_co2*(input$chamber_vol*input$atm_pres))/(R*(273.15+mean_temp)*input$chamber_area)

          if(var(data_subset$airt) == 0) {
            co2_flux = NA_real_
          } else {}

          results_string <- paste0(input$first_con_name,": slope = ", round(slope_ch4*3600, 2), " (ppm h<sup>-1</sup>)",
                                   ", flux = ", round(ch4_flux*3600, 2), " (\u03bcmol m<sup>-2</sup> h<sup>-1</sup>)",
                                   ", R<sup>2</sup> = ", round(r2_ch4, 2),
                                   "<br>",
                                   input$second_con_name,": slope = ", round(slope_co2*3600, 2), " (ppm h<sup>-1</sup>)",
                                   ", flux = ", round(co2_flux*3600, 2), " (\u03bcmol m<sup>-2</sup> h<sup>-1</sup>)",
                                   ", R<sup>2</sup> = ", round(r2_co2, 2))

          results <- data.frame("processing_date" = strftime(Sys.time(), "%Y-%m-%d %H:%M:%S", tz="GMT"),
                                "id" = as.character(input$sample_id),
                                "start" = strftime(ranges2$x[1], "%Y-%m-%d %H:%M:%S", tz="GMT"),
                                "end" = strftime(ranges2$x[2], "%Y-%m-%d %H:%M:%S", tz="GMT"),
                                "concentration1_name" = input$first_con_name,
                                "concentration1_slope_h1" = slope_ch4*3600,
                                "concentration1_intercept" = intercept_ch4,
                                "concentration1_R2" = r2_ch4,
                                "concentration2_name" = input$second_con_name,
                                "concentration2_slope_h1" = slope_co2*3600,
                                "concentration2_intercept" = intercept_co2,
                                "concentration2_R2" = r2_co2,
                                "temperature" = mean_temp,
                                "chamber_volume" = input$chamber_vol,
                                "chamber_area" = input$chamber_area,
                                "concentration1_flux_umol_m2_h1" = ch4_flux*3600,
                                "concentration2_flux_umol_m2_h1" = co2_flux*3600)
        }
      }else{
        data_subset <- data() %>%
          mutate(sec = cumsum(c(0, diff(as.numeric(datetime)))))

        results <- data.frame()
        results_string <- ""
      }
      return(list("df" = data_subset,
                  "results" = results,
                  "results_string" = results_string))

    })

    output$plot_zoom <- renderPlot({
      if(input$file_type == "Other") {
        validate(
          need(input$datetime_column, message = ""),
          need(input$concentration_values_ch4_column, message = ""))
      } else {}
      req(data_subset())
      zoom_data <- data_subset()
      req(input$first_con_name)
      zoom_plot_data <- zoom_data$df %>%
        filter(between(.data[[input$first_con_name]],input$ch4_range[1], input$ch4_range[2]),
               between(datetime, input$range[1], input$range[2]))

      if(nchar(input$second_con_name)>0){
        zoom_plot_data <- zoom_plot_data %>%
          filter(between(.data[[input$second_con_name]],input$co2_range[1],input$co2_range[2]))
      }

      par(mar = c(5,4,4,4) + 0.1)

      water_min = min(zoom_plot_data$water,na.rm=T)
      ch4_min = min(zoom_plot_data[[input$first_con_name]],na.rm=T)
      ch4_max = max(zoom_plot_data[[input$first_con_name]],na.rm=T)
      water_max = max(zoom_plot_data$water,na.rm=T)
      water_scaled = (ch4_max - ch4_min)*((zoom_plot_data$water-water_min)/(water_max - water_min))+ch4_min

      ggplot() +
        geom_point(data = zoom_plot_data, aes(.data$sec, .data[[input$first_con_name]], col = "CH4")) +
        geom_smooth(data = zoom_plot_data, aes(.data$sec, .data[[input$first_con_name]], col = "CH4"), method = "lm", se = F, formula = 'y ~ x') +
        geom_line(data = zoom_plot_data, aes(sec, water_scaled, col = "H2O")) +
        labs(x = "Time steps",
             y = input$first_con_name,
             col = "") +
        scale_color_manual(limits = c("CH4","H2O"),
                           labels = c(input$first_con_name,expression("H"[2]*"O")),
                           values = c("darkorange", "lightblue")) +
        theme_bw() +
        theme(legend.position = "bottom") -> p1

      if (nchar(input$second_con_name)==0) {
        print(p1)

      } else {

        co2_min = min(zoom_plot_data[[input$second_con_name]], na.rm=T)
        co2_max = max(zoom_plot_data[[input$second_con_name]], na.rm=T)

        co2_scaled = (ch4_max - ch4_min)*((zoom_plot_data[[input$second_con_name]]-co2_min)/(co2_max - co2_min))+ch4_min
        co2_labels = pretty(zoom_plot_data[[input$second_con_name]])
        co2_at = (ch4_max - ch4_min)*((co2_labels-co2_min)/(co2_max - co2_min))+ch4_min

        ch4_scaled = (co2_max - co2_min)*((zoom_plot_data[[input$first_con_name]]-ch4_min)/(ch4_max - ch4_min))+co2_min
        ch4_labels = pretty(zoom_plot_data[[input$first_con_name]])
        ch4_at = (co2_max - co2_min)*((ch4_labels-ch4_min)/(ch4_max - ch4_min))+co2_min

        if (!is.null(ranges2$x)){
        }

        lm_model_co2_scaled <- lm(co2_scaled~zoom_plot_data$sec, na.action=na.exclude)
        slope_co2_scaled <- coef(lm_model_co2_scaled)[2]
        intercept_co2_scaled <- coef(lm_model_co2_scaled)[1]

        p1 +
          geom_point(data = zoom_plot_data, aes(sec, co2_scaled, col = "CO2")) +
          geom_smooth(data = zoom_plot_data, aes(sec, co2_scaled, col = "CO2"), method = "lm", se = F,formula = 'y ~ x') +
          scale_color_manual(limits = c("CH4","CO2","H2O"),
                             labels = c(input$first_con_name,input$second_con_name,expression("H"[2]*"O")),
                             values = c("darkorange","forestgreen", "lightblue")) +
          scale_y_continuous(sec.axis = sec_axis(trans=~., name = input$second_con_name,
                                                 breaks = co2_at, labels = co2_labels))
      }
    })

    output$result_string <- renderUI({
      req(data_subset())
      if (is.null(ranges2$x)) return(NULL)
      HTML(data_subset()$results_string)
    })

    observe({
      brush <- input$plot_brush
      if (!is.null(brush)) {
        ranges2$x <- c(brush$xmin, brush$xmax)
        ranges2$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges2$x <- NULL
        ranges2$y <- NULL
      }
    })

    observeEvent(input$save,{
      data_out(rbind(data_out(), data_subset()$results) %>%
                 mutate(across(where(is.double), ~round(.x, 3))))
      output$results <- renderDT(
        data_out() %>%
          select(any_of(c("id", "start", "end","concentration1_name","concentration1_R2", "concentration1_flux_umol_m2_h1","concentration2_name","concentration2_R2", "concentration2_flux_umol_m2_h1"))),
        options = list(pageLength = 10, autoWidth = TRUE),
        rownames= FALSE)
    })

    observeEvent(input$delete_saved, {
      data_out(data.frame())

      output$results <- renderDT(
        data_out() %>%
          select(any_of(c("id", "start", "end","concentration1_name","concentration1_R2", "concentration1_flux_umol_m2_h1","concentration2_name","concentration2_R2", "concentration2_flux_umol_m2_h1"))),
        options = list(pageLength = 10, autoWidth = TRUE),
        rownames= FALSE)
    })


    # FluxSeparator ----
    ## Ebul ----
    ebul_prepared <- reactive({
      req(data())
      if (nchar(input$second_con_name) == 0) {
        data() %>%
          filter(between(datetime, input$range[1], input$range[2]),
                 between(.data[[input$first_con_name]], input$ch4_range[1], input$ch4_range[2])) %>%
          ungroup() %>%
          mutate(station = 1,
                 plot_number = floor((PumpCycle - min(PumpCycle)) / input$number_of_pumpcycles_in_plot) + 1,
                 sensor = 1)
      } else {
        data() %>%
          filter(between(datetime, input$range[1], input$range[2]),
                 between(.data[[input$first_con_name]], input$ch4_range[1], input$ch4_range[2]),
                 between(.data[[input$second_con_name]], input$co2_range[1], input$co2_range[2])) %>%
          ungroup() %>%
          mutate(station = 1,
                 plot_number = floor((PumpCycle - min(PumpCycle)) / input$number_of_pumpcycles_in_plot) + 1,
                 sensor = 1)
      }
    })

    observe({
      req(ebul_prepared())
      updateSliderInput(session, "plot_number_select",
                        min = min(ebul_prepared()$plot_number),
                        max = max(ebul_prepared()$plot_number), step = 1)
    })

    fluxsep_ebul <- reactive({
      req(ebul_prepared(), nzchar(input$concentration_values))
      ebul_prepared() %>%
        filter(plot_number == input$plot_number_select) %>%
        ebul_flux(concentration_values = input$concentration_values,
                  top_selection = input$top_selection,
                  IndexSpan = input$indexspan,
                  runvar_cutoff = input$runvar_cutoff,
                  concentration_diffusion_cutoff = input$concentration_diffusion_cutoff,
                  number_of_pumpcycles_in_plot = input$number_of_pumpcycles_in_plot,
                  smooth_data = input$smooth_data)
    }) |> bindCache(
      input$concentration_values, input$top_selection, input$indexspan,
      input$runvar_cutoff, input$concentration_diffusion_cutoff,
      input$number_of_pumpcycles_in_plot, input$smooth_data,
      input$range, input$ch4_range, input$co2_range, input$second_con_name,
      input$plot_number_select
    )

    observeEvent(input$concentration_values, {
      req(ebul_prepared())
      max_number_of_pump <- ebul_prepared() %>% ungroup() %>% distinct(PumpCycle) %>% nrow()
      updateSliderInput(session, "number_of_pumpcycles_in_plot", min = 1, max = max_number_of_pump, step = 1,
                        value = max_number_of_pump / 2)
    })

    output$plot2  <- renderPlot(fluxsep_ebul()$plot)

    output$ebul_table <- DT::renderDT(fluxsep_ebul()$bubbles %>%
                                        select(c(station:datetime_end, temp, concentration_per_hour = concentration_per_time)) %>%
                                        mutate(across(temp:concentration_per_hour, ~round(.x, 3)),
                                               datetime_start = strftime(datetime_start, format = '%Y-%m-%d %R'),
                                               datetime_end = strftime(datetime_end, format = '%Y-%m-%d %R')) %>%
                                        rename("concentration_per_hour (ppm h-1)" = concentration_per_hour,
                                               'airt (\u00b0C)' = temp),
                                      options = list(pageLength = 10, autoWidth = TRUE),
                                      rownames= FALSE)

    output$download_ebul <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_ebullitive_fluxdata.csv")
      },
      content = function(file) {
        fluxsep_ebul()$bubbles %>%
          rename("concentration_per_hour (ppm h-1)" = concentration_per_time,
                 'airt (\u00b0C)' = temp) %>%
          write.csv(file, row.names = FALSE)
      }
    )

    ## Diff ----
    diff_prepared <- reactive({
      req(data())
      if (nchar(input$second_con_name) == 0) {
        data() %>%
          filter(between(datetime, input$range[1], input$range[2]),
                 between(.data[[input$first_con_name]], input$ch4_range[1], input$ch4_range[2])) %>%
          ungroup() %>%
          mutate(station = "1",
                 plot_number = floor((PumpCycle - min(PumpCycle)) / input$number_of_pumpcycles_in_plot_diff) + 1,
                 sensor = "1")
      } else {
        data() %>%
          filter(between(datetime, input$range[1], input$range[2]),
                 between(.data[[input$first_con_name]], input$ch4_range[1], input$ch4_range[2]),
                 between(.data[[input$second_con_name]], input$co2_range[1], input$co2_range[2])) %>%
          ungroup() %>%
          mutate(station = "1",
                 plot_number = floor((PumpCycle - min(PumpCycle)) / input$number_of_pumpcycles_in_plot_diff) + 1,
                 sensor = "1")
      }
    })

    observe({
      req(diff_prepared(), nzchar(input$concentration_values_diff))
      df <- diff_prepared()
      cutoff_vals    <- df %>% pull(input$concentration_values_diff)
      max_antal_diff <- df %>% count(PumpCycle) %>% pull(n) %>% max()
      updateSliderInput(session, "plot_number_select_diff",
                        min = min(df$plot_number), max = max(df$plot_number), step = 1)
      updateSliderInput(session, "remove_observations_prior", min = 0, max = max_antal_diff, step = 1)
      updateSliderInput(session, "number_of_observations_used", min = 3, max = max_antal_diff, step = 1)
      updateSliderInput(session, "cutoff_start_value",
                        min = floor(min(cutoff_vals, na.rm = TRUE)),
                        max = ceiling(max(cutoff_vals, na.rm = TRUE)),
                        step = 1, value = ceiling(max(cutoff_vals, na.rm = TRUE)))
      updateSliderInput(session, "number_of_observations_required", min = 3, max = max_antal_diff, step = 1)
    })

    fluxsep_diff <- reactive({
      req(diff_prepared(), nzchar(input$concentration_values_diff))
      diff <- diff_prepared() %>%
        filter(plot_number == input$plot_number_select_diff) %>%
        diff_flux(concentration_values = input$concentration_values_diff,
                  IndexSpan = input$indexspan,
                  runvar_cutoff = input$runvar_cutoff,
                  remove_observations_prior = input$remove_observations_prior,
                  number_of_observations_used = input$number_of_observations_used,
                  cutoff_start_value = input$cutoff_start_value,
                  number_of_observations_required = input$number_of_observations_required,
                  number_of_pumpcycles_in_plot = input$number_of_pumpcycles_in_plot_diff,
                  smooth_data = input$smooth_data_diff,
                  look_for_bubbles = input$look_for_bubbles,
                  Hutchinson_Mosier_correction = input$hmr_correction,
                  volume = input$volume_diff,
                  area = input$area_diff)
      validate(
        need(diff[1] != "ERROR! NO DATA",
             message = "Couldn't find any diffusive fluxes. Try lowering the runvar_cutoff or turning off looking for bubbles")
      )
      diff
    }) |> bindCache(
      input$concentration_values_diff, input$indexspan, input$runvar_cutoff,
      input$remove_observations_prior, input$number_of_observations_used,
      input$cutoff_start_value, input$number_of_observations_required,
      input$number_of_pumpcycles_in_plot_diff, input$smooth_data_diff,
      input$look_for_bubbles, input$hmr_correction, input$volume_diff, input$area_diff,
      input$range, input$ch4_range, input$co2_range, input$second_con_name,
      input$plot_number_select_diff
    )

    observeEvent(input$concentration_values_diff, {
      req(diff_prepared())
      max_number_of_pump_diff <- diff_prepared() %>% ungroup() %>% distinct(PumpCycle) %>% nrow()
      updateSliderInput(session, "number_of_pumpcycles_in_plot_diff", min = 1, max = max_number_of_pump_diff, step = 1,
                        value = max_number_of_pump_diff / 2)
    })

    output$plot_diff  <- renderPlot(fluxsep_diff()$plot)

    output$diff_table <- renderDT(fluxsep_diff()$diff %>%
                                    mutate(across(slope_concentration_hr:temp, ~round(.x,3)),
                                           datetime_start = strftime(datetime_start, format = '%Y-%m-%d %R'),
                                           datetime_end = strftime(datetime_end, format = '%Y-%m-%d %R')) %>%
                                    rename("slope_concentration_hr (ppm h-1)" = slope_concentration_hr,
                                           'airt (\u00b0C)' = temp),
                                  options = list(pageLength = 10, autoWidth = TRUE),
                                  rownames= FALSE)

    output$download_diff <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_diffusive_fluxdata.csv")
      },
      content = function(file) {
        fluxsep_diff()$diff %>%
          rename("slope_concentration_hr (ppm h-1)" = slope_concentration_hr,
                 'airt (\u00b0C)' = temp) %>%
          write.csv(file, row.names = FALSE)
      }
    )

    # Downloads ----

    output$download <- downloadHandler(

      filename = function() {
        paste0(Sys.Date(), "_fluxdata.csv")
      },

      content = function(file) {
        write.csv(data_out(), file, row.names = FALSE)
      })

    output$ch4_download <- downloadHandler(

      filename = function() {
        paste0(Sys.Date(), "_data.csv")
      },

      content = function(file) {
        data() %>%
          rename('rh (%)' = rh, 'airt (\u00b0C)' = airt, "H2O (ppm)" = water)  %>%
          write.csv(file, row.names = FALSE)
      })

  }

  shinyApp(ui, server, ...)
}
