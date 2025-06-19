rm(list = ls(all = TRUE))
library(ncdf4)
library(raster)
library(udunits2)
library(lpjmlkit)

# Directory containing NetCDF source data from ISIMIP.
isimip_source_dir <- "SOURCE"
isimip_source_unit <- "g/m2/month"

# Directory where generated input files for LPJmL are saved.
output_dir <- "TARGET"

# A grid file in LPJmL format. Deposition inputs are generated for this grid.
gridfilename <- "GRIDFILENAME"

# CLM version created. Version 1 & 2 use data type SHORT, version 3 uses FLOAT.
# FLOAT allows higher precision and is therefore preferable.
depo_version <- 3
# Scaling factor used for deposition files in case of lu_version < 3
depo_scalar <- 0.001
depo_order <- 1
force_overwrite <- FALSE # if set to FALSE, will skip conversion if output file exists
lpj_unit <- "g/m2/day"

# Time series in ISIMIP source data is always split into two files.
scenario_years <- list(
  "1901soc" = list(1850:1900, 1901:2021),
  "2015soc" = list(1850:1900, 1901:2021),
  histsoc = list(1850:1900, 1901:2021),
  "1850soc" = list(1850:1900, 1901:2021)
)
scenario_gcm_specific <- list(
  "1901soc" = FALSE,
  "2015soc" = FALSE,
  histsoc = FALSE,
  "1850soc" = FALSE
)

output_years <- list(
  "1901soc" = 1850:2021,
  "2015soc" = 1850:2021,
  histsoc = 1850:2021,
  "1850soc" = 1850:2021
)

output_scenarios <- list(
  "1901soc" = "1901soc",
  "2015soc" = "2015soc",
  histsoc = "histsoc",
  "1850soc" = "1850soc"
)

output_gcm_specific <- list(
  "1901soc" = FALSE,
  "2015soc" = FALSE,
  histsoc = FALSE,
  "1850soc" = FALSE
)

# GCMs not relevant in ISIMIP3a
GCMs <- c(
  "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"
)
GCMs_secondary <- c(
  "CanESM5", "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "MIROC6"
)

ndayspermonth <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
ndaysperyear <- 365

if (file.exists(gridfilename)) {
  gridheader <- lpjmlkit::read_header(gridfilename)
  griddata <- lpjmlkit::read_grid(gridfilename)$data
  gridarea <- lpjmlkit:: calc_cellarea(
    griddata[, "lat"],
    gridheader$header["cellsize_lon"],
    gridheader$header["cellsize_lat"],
    return_unit = "m2"
  )
} else {
  stop(paste("Error: cannot find grid file", gridfilename))
}

# Determine required unit conversion.
mon2day <- year2day <- FALSE
if (grepl("/mon|mon[[:alpha:]]*-1", isimip_source_unit)) {
  timeunit <- regmatches(
    isimip_source_unit,
    regexpr("(/)*mon[[:alpha:]]*(-1)*", isimip_source_unit)
  )
  mon2day <- TRUE
  isimip_conv_unit <- sub(timeunit, "/day", isimip_source_unit)
  cat("Converting source data from time unit", timeunit, "to /day\n")
} else if (grepl("/yr|/year|yr-1|year-1", isimip_source_unit)) {
  timeunit <- regmatches(
    isimip_source_unit,
    regexpr("/yr|/year|yr-1|year-1", isimip_source_unit)
  )
  year2day <- TRUE
  isimip_conv_unit <- sub(timeunit, "/day", isimip_source_unit)
  cat("Converting source data from time unit", timeunit, "to /day\n")
} else if (grepl("day", isimip_source_unit)) {
  isimip_conv_unit <- isimip_source_unit
} else {
  stop("Unknown source unit")
}
if (!udunits2::ud.are.convertible(isimip_conv_unit, lpj_unit)) {
  stop("Unit conversion error")
}
for (output in names(output_scenarios)) {
  print(paste("**** Process", output, "****"))
  for (gcm in switch(
    as.character(output_gcm_specific[[output]]),
    "TRUE" = GCMs,
    "FALSE" = ""
  )) {
    if (gcm != "") {
      print(paste("++++ Process", gcm, "++++"))
    }
    fileyears <- range(output_years[[output]])
    inputfiles_no3 <- list()
    inputfiles_nh4 <- list()
    cell_index_no3 <- cell_index_nh4 <- list()
    avail_years_no3 <- avail_years_nh4 <- integer()
    for (scen in output_scenarios[[output]]) {
      if (is.list(scenario_years[[scen]])) {
        inputname_no3 <- file.path(
          isimip_source_dir, scen,
          paste0(
            "ndep-noy_", scen, "_monthly_",
            sapply(scenario_years[[scen]], min), "_",
            sapply(scenario_years[[scen]], max),
            ".nc"
          )
        )
        inputname_nh4 <- file.path(
          isimip_source_dir, scen,
          paste0(
            "ndep-nhx_", scen, "_monthly_",
            sapply(scenario_years[[scen]], min), "_",
            sapply(scenario_years[[scen]], max),
            ".nc"
          )
        )
      } else {
        inputname_no3 <- file.path(
          isimip_source_dir, scen,
          paste0(
            "ndep-noy_", scen, "_monthly_",
            min(scenario_years[[scen]]), "_",
            max(scenario_years[[scen]]),
            ".nc"
          )
        )
        inputname_nh4 <- file.path(
          isimip_source_dir, scen,
          paste0(
            "ndep-nhx_", scen, "_monthly_",
            min(scenario_years[[scen]]), "_",
            max(scenario_years[[scen]]),
            ".nc"
          )
        )
      }
      for (f in seq_along(inputname_no3)) {
        if (file.exists(inputname_no3[f])) {
          inputfiles_no3[[scen]][[f]] <- ncdf4::nc_open(inputname_no3[f])
          tmpraster <- raster::raster(inputname_no3[f], band = 1)
          cell_index_no3[[scen]][[f]] <- raster::cellFromXY(tmpraster, griddata)
          if (is.list(scenario_years[[scen]])) {
            avail_years_no3 <- union(avail_years_no3, scenario_years[[scen]][[f]])
          } else {
            avail_years_no3 <- union(avail_years_no3, scenario_years[[scen]])
          }
        } else {
          print(
            paste0(
              "Warning: expected file ", inputname_no3[f],
              " not found. Removing scenario_years[[", scen, "]][[", f, "]]"
            )
          )
          inputfiles_no3[[scen]][[f]] <- FALSE
          # Remove available years because file does not exist.
          if (is.list(scenario_years[[scen]])) {
            scenario_years[[scen]][[f]] <- NULL
          } else {
            scenario_years[[scen]] <- NULL
          }
        }
        if (file.exists(inputname_nh4[f])) {
          inputfiles_nh4[[scen]][[f]] <- ncdf4::nc_open(inputname_nh4[f])
          tmpraster <- raster::raster(inputname_nh4[f], band = 1)
          cell_index_nh4[[scen]][[f]] <- raster::cellFromXY(tmpraster, griddata)
          if (is.list(scenario_years[[scen]])) {
            avail_years_nh4 <- union(avail_years_nh4, scenario_years[[scen]][[f]])
          } else {
            avail_years_nh4 <- union(avail_years_nh4, scenario_years[[scen]])
          }
        } else {
          print(
            paste0(
              "Warning: expected file ", inputname_nh4[f],
              " not found. Removing scenario_years[[", scen, "]][[", f, "]]"
            )
          )
          inputfiles_nh4[[scen]][[f]] <- FALSE
          # Remove available years because file does not exist.
          if (is.list(scenario_years[[scen]])) {
            scenario_years[[scen]][[f]] <- NULL
          } else {
            scenario_years[[scen]] <- NULL
          }
        }
      }
    }
    if (!all(seq(min(fileyears), max(fileyears)) %in% avail_years_no3)) {
      possible <- intersect(seq(min(fileyears), max(fileyears)), avail_years_no3)
      if (length(possible) > 1 &&
        length(possible) != max(possible) - min(possible) + 1
      ) {
        # avail_years has breaks
        for (i in seq(2, length(possible))) {
          if (possible[i] - possible[i - 1] != 1) {
            possible <- possible[seq(1, i - 1)]
            break
          }
        }
      }
      print(
        paste0(
          "Warning: available NO3 source files do not cover the full requested ",
          "output_years ",
          paste(range(output_years[[output]]), collapse = "-"),
          ". Shortening output file to ",
          paste(range(possible), collapse = "-")
        )
      )
      fileyears <- range(possible)
    }
    if (!all(seq(min(fileyears), max(fileyears)) %in% avail_years_nh4)) {
      possible <- intersect(seq(min(fileyears), max(fileyears)), avail_years_nh4)
      if (length(possible) > 1 &&
        length(possible) != max(possible) - min(possible) + 1
      ) {
        # avail_years has breaks
        for (i in seq(2, length(possible))) {
          if (possible[i] - possible[i - 1] != 1) {
            possible <- possible[seq(1, i - 1)]
            break
          }
        }
      }
      print(
        paste0(
          "Warning: available NH4 source files do not cover the full requested ",
          "output_years ",
          paste(range(output_years[[output]]), collapse = "-"),
          ". Shortening output file to ",
          paste(range(possible), collapse = "-")
        )
      )
      fileyears <- range(possible)
    }
    # LPJmL filename
    filename_no3 <- file.path(
      output_dir,
      paste0(
        "no3-deposition_isimip3a_", output, "_",
        ifelse(output_gcm_specific[[output]], paste0(gcm, "_"), ""),
        fileyears[1], "-", fileyears[2], ".clm"
      )
    )
    filename_nh4 <- file.path(
      output_dir,
      paste0(
        "nh4-deposition_isimip3a_", output, "_",
        ifelse(output_gcm_specific[[output]], paste0(gcm, "_"), ""),
        fileyears[1], "-", fileyears[2], ".clm"
      )
    )
    depoheader <- lpjmlkit::create_header(
      name = "LPJNDEP",
      version = depo_version,
      order = depo_order,
      firstyear = min(fileyears),
      nyear = max(fileyears) - min(fileyears) + 1,
      firstcell = 0,
      ncell = gridheader$header["ncell"],
      nbands = 12,
      cellsize_lon = gridheader$header["cellsize_lon"],
      scalar = ifelse(depo_version < 3, depo_scalar, 1),
      cellsize_lat = gridheader$header["cellsize_lat"],
      datatype = ifelse(depo_version < 3, 1, 3),
      endian = .Platform$endian
    )
    fe <- file.exists(c(filename_nh4, filename_no3))
    names(fe) <- c(filename_nh4, filename_no3)
    if (any(fe) && !force_overwrite) {
      message(
        "Skipping ", output, " because output file(s) ",
        toString(names(fe)[which(fe)]), " exist(s)"
      )
      expected_size <- lpjmlkit::get_headersize(depoheader) + 
        prod(depoheader$header[c("nbands", "ncell", "nyear")]) *
        lpjmlkit::get_datatype(depoheader)$size
      fs <- file.size(c(filename_nh4, filename_no3))
      names(fs) <- c(filename_nh4, filename_no3)
      if (any(fs != expected_size)) {
        warning(
          "file(s) ", toString(names(fs)[which(fs != expected_size)]),
          " has/have unexpected size. Expected: ",
          expected_size, ", detected: ",
          toString(fs[which(fs != expected_size)]),
          ". Please delete existing file(s) ",
          "set force_overwrite to TRUE.",
          call. = FALSE, immediate. = TRUE
        )
      }
      sapply(inputfiles_nh4[[scen]], ncdf4::nc_close)
      sapply(inputfiles_no3[[scen]], ncdf4::nc_close)
      next
    }
    lpjmlkit::write_header(filename_nh4, depoheader, force_overwrite)
    lpjmlkit::write_header(filename_no3, depoheader, force_overwrite)
    outfile_nh4 <- file(filename_nh4, "ab")
    outfile_no3 <- file(filename_no3, "ab")
    cat("Converting", depoheader$header["nyear"], "years\n")
    for (y in seq(min(fileyears), max(fileyears))) {
      if (y %% 10 == 0)
        print(y)
      # find which file to use for this year
      for (check_scen in output_scenarios[[output]]) {
        if (is.list(scenario_years[[check_scen]])) {
          for (f in seq_along(scenario_years[[check_scen]])) {
            if (y %in% scenario_years[[check_scen]][[f]]) {
              # pointer to NetCDF file
              nh4_nc <- inputfiles_nh4[[check_scen]][[f]]
              no3_nc <- inputfiles_no3[[check_scen]][[f]]
              nc2grid_nh4 <- cell_index_nh4[[check_scen]][[f]]
              nc2grid_no3 <- cell_index_no3[[check_scen]][[f]]
              # find year in file
              time_index <- (y - min(scenario_years[[check_scen]][[f]])) * 12 + 1
              break
            }
          }
        } else {
          if (y %in% scenario_years[[check_scen]]) {
            # pointer to NetCDF file
            nh4_nc <- inputfiles_nh4[[check_scen]][[1]]
            no3_nc <- inputfiles_no3[[check_scen]][[1]]
            nc2grid_nh4 <- cell_index_nh4[[check_scen]][[1]]
            nc2grid_no3 <- cell_index_no3[[check_scen]][[1]]
            # find year in file
            time_index <- (y - min(scenario_years[[check_scen]])) * 12 + 1
            break
          }
        }
      }
      if (!exists("nh4_nc")) {
        stop("Unexpected error")
      }
      nh4_flip <- (nh4_nc$dim$lat$vals[1] < nh4_nc$dim$lat$vals[2])
      no3_flip <- (no3_nc$dim$lat$vals[1] < no3_nc$dim$lat$vals[2])
      # array for LPJ data
      nh4_year <- no3_year <- array(
        0,
        dim = c(depoheader$header["nbands"], depoheader$header["ncell"])
      )
      nh4_year_nc <- ncdf4::ncvar_get(
        nh4_nc, "nhx", start = c(1, 1, time_index), count = c(-1, -1, 12)
      )
      if (nh4_flip) {
        nh4_year_nc <- nh4_year_nc[, seq(nh4_nc$dim$lat$len, 1), ]
      }
      no3_year_nc <- ncdf4::ncvar_get(
        no3_nc, "noy", start = c(1, 1, time_index), count = c(-1, -1, 12)
      )
      if (no3_flip) {
        no3_year_nc <- no3_year_nc[, seq(no3_nc$dim$lat$len, 1), ]
      }
      dim(nh4_year_nc) <- c(prod(dim(nh4_year_nc)[c(1, 2)]), dim(nh4_year_nc)[3])
      dim(no3_year_nc) <- c(prod(dim(no3_year_nc)[c(1, 2)]), dim(no3_year_nc)[3])
      nh4_year <- t(nh4_year_nc[nc2grid_nh4, ])
      no3_year <- t(no3_year_nc[nc2grid_no3, ])
      if (mon2day) {
        nh4_year <- nh4_year / ndayspermonth
        no3_year <- no3_year / ndayspermonth
      } else if (year2day) {
        nh4_year <- nh4_year / ndaysperyear
        no3_year <- no3_year / ndaysperyear
      }
      nh4_year <- udunits2::ud.convert(nh4_year, isimip_conv_unit, lpj_unit)
      no3_year <- udunits2::ud.convert(no3_year, isimip_conv_unit, lpj_unit)
      if (depo_order == 4) {
        nh4_year <- t(nh4_year)
        no3_year <- t(no3_year)
      }
      if (typeof(lpjmlkit::get_datatype(depoheader)$type) == "integer") {
        max_range <- 2^(lpjmlkit::get_datatype(depoheader)$size * 8 - 1)
        if (any(round(nh4_year / depoheader$header["scalar"]) > max_range) ||
          any(round(no3_year / depoheader$header["scalar"]) > max_range)
        ) {
          excess <- c(
            length(which(round(nh4_year / depoheader$header["scalar"]) > max_range)),
            length(which(round(no3_year / depoheader$header["scalar"]) > max_range))
          )
          names(excess) <- c(nh4_nc$filename, no3_nc$filename)
          warning(
            "Year ", y, ": ",
            toString(excess[which(excess > 0)]), " cell(s) in file(s) ",
            toString(names(which(excess > 0))),
            " above maximum value range of data type. Limiting to valid range. ",
            "Consider changing scalar or switching to float data type.",
            call. = FALSE, immediate. = TRUE
          )
          nh4_year <- pmin(
            nh4_year / depoheader$header["scalar"] , max_range
          ) * depoheader$header["scalar"]
          no3_year <- pmin(
            no3_year / depoheader$header["scalar"] , max_range
          ) * depoheader$header["scalar"]
        }
        nh4_year <- as.integer(round(nh4_year / depoheader$header["scalar"]))
        no3_year <- as.integer(round(no3_year / depoheader$header["scalar"]))
      } else {
        nh4_year <- as.double(nh4_year / depoheader$header["scalar"])
        no3_year <- as.double(no3_year / depoheader$header["scalar"])
      }
      writeBin(
        nh4_year,
        outfile_nh4,
        size = lpjmlkit::get_datatype(depoheader)$size,
        endian = depoheader$endian
      )
      writeBin(
        no3_year,
        outfile_no3,
        size = lpjmlkit::get_datatype(depoheader)$size,
        endian = depoheader$endian
      )
      rm(nh4_nc, no3_nc)
    }
    close(outfile_nh4)
    close(outfile_no3)
    for (scen in names(inputfiles_nh4))
      sapply(inputfiles_nh4[[scen]], ncdf4::nc_close)
    for (scen in names(inputfiles_no3))
      sapply(inputfiles_no3[[scen]], ncdf4::nc_close)
    cat(
      "Finished writing", depoheader$header["nyear"],
      "years to", toString(c(filename_nh4, filename_no3)), "\n"
    )
  }
}
