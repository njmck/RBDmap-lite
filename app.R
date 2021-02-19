#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Run the following in the console if not all packages are installed/updated:
# install.packages(c("shiny", "seqinr", "tidyverse", "openxlsx", "plyr", "RColorBrewer", "rlist"))

library(shiny)
library(seqinr)
library(tidyverse)
library(openxlsx)
library(plyr)
library(RColorBrewer)
library(rlist)

# Define UI to allow users to upload their files.
ui <- fluidPage(
		title = "R Shiny - RBDmap Lite",
		titlePanel(
			h1("RBDmap Lite", align = "center")
			),
  		
		sidebarLayout(
			sidebarPanel(
				fileInput(inputId = "fasta_file", multiple = FALSE, label = "Select fasta file:", buttonLabel = "Browse..."),
				fileInput(inputId = "upload_files", multiple = TRUE, label = "Select RBDmap files (tab-separated text files, multiple okay):", buttonLabel = "Browse..."),
				actionButton(inputId = "analyse_button", label = "Analyse!")
				),		
			mainPanel(
				uiOutput("prot_select"),
				tableOutput("filedf"),
				textOutput(outputId = "filedf_text"),
				plotOutput(outputId = "repeats_plot"),
				htmlOutput(outputId = "fasta_color_legend"),
				htmlOutput(outputId = "fasta_color"),
				tableOutput("total_df")
		  		)
		)
)

# Server
server <- function(input, output) {
    # When the user presses the button, the following code runs:
	observeEvent(input$analyse_button, {
    
		# R code to process user uploaded files:    
    	protein_fasta <- read.fasta(input$fasta_file[1,'datapath'], seqtype = "AA")
    	target_proteins <- getName(protein_fasta)
		
    	# User can select by drop-down menu which protein within the fasta to analyse:
    	output$prot_select <- renderUI({
        	selectInput(inputId = "fasta_select",
            	choices = sort(target_proteins),
                label = "Select:")
    	})
    	
    	filelist <- input$upload_files[,'name']
    	RBDmap_tibble <- NULL
    	
    	for (rows in 1:length(input$upload_files[,"datapath"])) {
        	tibble <- read_tsv(input$upload_files[rows,"datapath"])
        	tibble <- add_column(tibble, Filename = input$upload_files[rows,"name"])
        	RBDmap_tibble <- rbind.fill(RBDmap_tibble, tibble)
    	}
    	
    	# Create new columns for sample type, and fill ----
    	
    	RBDmap_tibble <- add_column(RBDmap_tibble, Treatment = "")
    	RBDmap_tibble <- add_column(RBDmap_tibble, Sample = "")
    	RBDmap_tibble <- add_column(RBDmap_tibble, Protease = "")
    	RBDmap_tibble <- add_column(RBDmap_tibble, Repeat_No = "")
    	
    	for (row in 1:nrow(RBDmap_tibble)) {
        	# NoXL or UV:
        	if (str_detect(RBDmap_tibble[row, "Filename"], regex("NoXL", ignore_case = TRUE)) == TRUE) {
        	    RBDmap_tibble[row, "Treatment"] <- "NoXL"
        	} else if (str_detect(RBDmap_tibble[row, "Filename"], regex("UV", ignore_case = TRUE)) == TRUE) {
        	    RBDmap_tibble[row, "Treatment"] <- "UV"
        	} else {
        	    #print("ERROR: NO 'TREATMENT' FOUND IN FILENAME")
        	}
        	# Input or Eluate:
        	if (str_detect(RBDmap_tibble[row, "Filename"], regex("Input", ignore_case = TRUE)) == TRUE) {
        	    RBDmap_tibble[row, "Sample"] <- "Input"
        	} else if (str_detect(RBDmap_tibble[row, "Filename"], regex("Eluate", ignore_case = TRUE)) == TRUE) {
        	    RBDmap_tibble[row, "Sample"] <- "Eluate"
        	} else {
        	    #print("ERROR: NO 'SAMPLE' FOUND IN FILENAME")
        	}
        	# ArgC or LysC:
        	if (str_detect(RBDmap_tibble[row, "Filename"], regex("ArgC", ignore_case = TRUE)) == TRUE) {
        	    RBDmap_tibble[row, "Protease"] <- "ArgC"
        	} else if (str_detect(RBDmap_tibble[row, "Filename"], regex("LysC", ignore_case = TRUE)) == TRUE) {
        	    RBDmap_tibble[row, "Protease"] <- "LysC"
        	} else {
        	    #print("ERROR: NO 'PROTEASE' FOUND IN FILENAME")
        	}
        	# Repeat Number column:
        	repeat_no <- as.numeric(paste(unlist(str_match_all(RBDmap_tibble[row, "Filename"], "[[:digit:]]")), collapse = ""))
        	if (repeat_no != "") {
        	    RBDmap_tibble[row, "Repeat_No"] <- repeat_no
        	} else {
        	    #print("ERROR: NO 'REPEAT_NO' FOUND IN FILENAME")
        	}
		}
    	
    	# Create separate "Input" and "Eluate" tibbles ----
    	
    	input_tibble <- NULL
    	
    	for (row in 1:nrow(RBDmap_tibble)) {
        	if (RBDmap_tibble[row, "Sample"] == "Input") {
            	input_tibble <- bind_rows(input_tibble, RBDmap_tibble[row,])
        	}
		}
    	
    	eluate_tibble <- NULL
    	
    	for (row in 1:nrow(RBDmap_tibble)) {
        	if (RBDmap_tibble[row, "Sample"] == "Eluate") {
            	eluate_tibble <- bind_rows(eluate_tibble, RBDmap_tibble[row,])
        	}
    	}
    	
    	# Delete peptides from eluate tibble that do not contain anything from target_proteins ----
    	
    	delete_list_1 <- c()
    	
    	for (row in 1:nrow(eluate_tibble)) {
        	if (eluate_tibble[row,"Proteins"] %in% target_proteins) {
            	#print("MATCH FOR TARGET")
        	} else {
            	#print("NOT MATCH")
            	delete_list_1 <- c(delete_list_1, row)
        	}
    	}
    	
    	eluate_tibble <- eluate_tibble[-delete_list_1,]
    	
    	# Delete peptides below the score threshold ----
    	
    	andromeda_threshold <- 20 # Set your own Andromeda score threshold here.
    	
    	delete_list_2 <- c()
    	
    	for (row in 1:nrow(eluate_tibble)) {
        	if (max(eluate_tibble[row,"Score"], andromeda_threshold) == eluate_tibble[row,"Score"]) {
            	#print(paste(toString(eluate_tibble[row,"Score"]),"IS HIGHER THAN", toString(andromeda_threshold)))
        	} else {
            	#print(paste(toString(eluate_tibble[row,"Score"]),"IS LOWER THAN", toString(andromeda_threshold)))
            	delete_list_2 <- c(delete_list_2, row)
        	}
    	}
    	
    	eluate_tibble <- eluate_tibble[-delete_list_2,]
    	
    	# Delete peptides that don't start and end with K/R ----
    	
    	delete_list_3 <- c()
    	
    	for (row in 1:nrow(eluate_tibble)) {
        	if (eluate_tibble[row,4] == eluate_tibble[row,8]) {
            	#print(paste(toString(eluate_tibble[row,4])," is the same as ",toString(eluate_tibble[row,8])))
            	delete_list_3 <- c(delete_list_3, row)
        	} else {
            	#print(paste(toString(eluate_tibble[row,4])," is NOT the same as ",toString(eluate_tibble[row,8])))
        	}
    	}
    	
    		if (is.null(delete_list_3) == TRUE) {
        		#print("NO ROWS TO DELETE")
    		} else {
    	    	eluate_tibble <- eluate_tibble[-delete_list_3,]
    	    	#print("DELETED ROWS")
    		}
    	
    	# Getting the cross-linked peptides ----
    
		eluate_tibble <- add_column(eluate_tibble, Crosslinked_Peptide = "")
		eluate_tibble <- add_column(eluate_tibble, XL_pep_start = "")
		eluate_tibble <- add_column(eluate_tibble, XL_pep_end = "")
    
    	# ArgC-digested peptides, with K before
    
    	for (row in 1:nrow(eluate_tibble)) {
    	    if ((length(grep("[[:alpha:]]", eluate_tibble[row,"N-term cleavage window"])) == 0)
    	        & (eluate_tibble[row,"Protease"] == "ArgC")
    	        & (eluate_tibble[row,"Amino acid before"] == "K")) {
    	        #print(paste("The sequence", toString(eluate_tibble[row,"Sequence"]), "is INCOMPATIBLE"))
    	    } else if ((eluate_tibble[row,"Protease"] == "ArgC") & (eluate_tibble[row,"Amino acid before"] == "K")) {
    	        protein_of_interest <- toString(eluate_tibble[row,"Proteins"])
    	        last_R <- as.numeric((max(grep("R", protein_fasta[[protein_of_interest]][1:(as.numeric(eluate_tibble[row,"Start position"]) - 2)])) + 1))
    	        start_pos <- as.numeric((eluate_tibble[row,"Start position"]) - 1)
    	        crosslinked_peptide <- paste(protein_fasta[[protein_of_interest]][last_R:start_pos],collapse = '')
    	        #print(paste("For the sequence:", toString(eluate_tibble[row,"Sequence"])))
    	        #print(crosslinked_peptide)
    	        eluate_tibble[row,"XL_pep_start"] <- min(last_R, start_pos)
    	        eluate_tibble[row,"XL_pep_end"] <- max(last_R, start_pos)
    	        eluate_tibble[row,"Crosslinked_Peptide"] <- crosslinked_peptide
    	    }
    	}
    
    	# ArgC-digested peptides, with R before
    
    	for (row in 1:nrow(eluate_tibble)) {
        	if ((length(grep("[[:alpha:]]", eluate_tibble[row,"C-term cleavage window"])) == 0)
            	& (eluate_tibble[row,"Protease"] == "ArgC")
            	& (eluate_tibble[row,"Amino acid before"] == "R")) {
            	#print(paste("The sequence", toString(eluate_tibble[row,"Sequence"]), "is INCOMPATIBLE"))
        	} else if ((eluate_tibble[row,"Protease"] == "ArgC")
        	           & (eluate_tibble[row,"Amino acid before"] == "R")) {
        	    protein_of_interest <- toString(eluate_tibble[row,"Proteins"])
        	    next_R <- as.numeric(min(grep("R", protein_fasta[[protein_of_interest]][as.numeric(eluate_tibble[row,"End position"]):length(protein_fasta[[protein_of_interest]])])) + (as.numeric(eluate_tibble[row,"End position"]) - 1))
        	    start_pos <- as.numeric((eluate_tibble[row,"End position"]) + 1)
        	    crosslinked_peptide <- paste(protein_fasta[[protein_of_interest]][start_pos:next_R],collapse = '')
        	    #print(paste("For the sequence:", toString(eluate_tibble[row,"Sequence"])))
        	    #print(crosslinked_peptide)
        	    eluate_tibble[row,"XL_pep_start"] <- min(next_R, start_pos)
        	    eluate_tibble[row,"XL_pep_end"] <- max(next_R, start_pos)
        	    eluate_tibble[row,"Crosslinked_Peptide"] <- crosslinked_peptide
        	}
    	}
    
    	# LysC-digested peptides, with R before
    	
    	for (row in 1:nrow(eluate_tibble)) {
    	    if ((length(grep("[[:alpha:]]", eluate_tibble[row,"N-term cleavage window"])) == 0)
    	        & (eluate_tibble[row,"Protease"] == "LysC")
    	        & (eluate_tibble[row,"Amino acid before"] == "R")) {
    	        #print(paste("The sequence", toString(eluate_tibble[row,"Sequence"]), "is INCOMPATIBLE"))
    	    } else if ((eluate_tibble[row,"Protease"] == "LysC") & (eluate_tibble[row,"Amino acid before"] == "R")) {
    	        protein_of_interest <- toString(eluate_tibble[row,"Proteins"])
    	        last_R <- as.numeric((max(grep("K", protein_fasta[[protein_of_interest]][1:(as.numeric(eluate_tibble[row,"Start position"]) - 2)])) + 1))
    	        start_pos <- as.numeric((eluate_tibble[row,"Start position"]) - 1)
    	        crosslinked_peptide <- paste(protein_fasta[[protein_of_interest]][last_R:start_pos],collapse = '')
    	        #print(paste("For the sequence:", toString(eluate_tibble[row,"Sequence"])))
    	        #print(crosslinked_peptide)
    	        eluate_tibble[row,"XL_pep_start"] <- min(last_R, start_pos)
    	        eluate_tibble[row,"XL_pep_end"] <- max(last_R, start_pos)
    	        eluate_tibble[row,"Crosslinked_Peptide"] <- crosslinked_peptide
    	    }
    	}
    
    	# LysC-digested peptides, with K before
    	
    	for (row in 1:nrow(eluate_tibble)) {
    	    if ((length(grep("[[:alpha:]]", eluate_tibble[row,"C-term cleavage window"])) == 0)
    	        & (eluate_tibble[row,"Protease"] == "LysC")
    	        & (eluate_tibble[row,"Amino acid before"] == "K")) {
    	        #print(paste("The sequence", toString(eluate_tibble[row,"Sequence"]), "is INCOMPATIBLE"))
    	    } else if ((eluate_tibble[row,"Protease"] == "LysC")
    	               & (eluate_tibble[row,"Amino acid before"] == "K")) {
    	        protein_of_interest <- toString(eluate_tibble[row,"Proteins"])
    	        next_R <- as.numeric(min(grep("K", protein_fasta[[protein_of_interest]][as.numeric(eluate_tibble[row,"End position"]):length(protein_fasta[[protein_of_interest]])])) + (as.numeric(eluate_tibble[row,"End position"]) - 1))
    	        start_pos <- as.numeric((eluate_tibble[row,"End position"]) + 1)
    	        crosslinked_peptide <- paste(protein_fasta[[protein_of_interest]][start_pos:next_R],collapse = '')
    	        #print(paste("For the sequence:", toString(eluate_tibble[row,"Sequence"])))
    	        #print(crosslinked_peptide)
    	        eluate_tibble[row,"XL_pep_start"] <- min(next_R, start_pos)
    	        eluate_tibble[row,"XL_pep_end"] <- max(next_R, start_pos)
    	        eluate_tibble[row,"Crosslinked_Peptide"] <- crosslinked_peptide
    	    }
    	}
    
    	# Get rid of all peptides for which no cross-linked peptide was found
    	
    	eluate_tibble <- filter(eluate_tibble, Crosslinked_Peptide != "")
    
    	# Make hit graphs for each protein ----
    	
    	table_list <- c()
    	plot_list <- list()
    	
    	for (prot_seq in 1:length(target_proteins)) {
        	graph_tibble <- NULL
        	graph_tibble <- tibble(residue = protein_fasta[[prot_seq]], res_no = 1:length(protein_fasta[[prot_seq]]), freq = 0)
        	# Loop in a loop:
        	for (row in 1:nrow(eluate_tibble)) {
        	    if (eluate_tibble[row, "Proteins"] == target_proteins[prot_seq]) {
        	        graph_tibble[as.numeric(eluate_tibble[row,"XL_pep_start"]):as.numeric(eluate_tibble[row,"XL_pep_end"]), "freq"] <- graph_tibble[as.numeric(eluate_tibble[row,"XL_pep_start"]):as.numeric(eluate_tibble[row,"XL_pep_end"]), "freq"] + 1
        	    }
        	}
        	assign(toString(make.names(paste(toString(target_proteins[prot_seq]), "_graph_tibble", sep = ""))), graph_tibble)
        	table_list <- c(table_list, toString(make.names(paste(toString(target_proteins[prot_seq]), "_graph_tibble", sep = ""))))
        	plot <- ggplot(graph_tibble, aes(res_no, freq)) + geom_step() +
        			labs(title = paste("Number of RBDmap Hits per Residue", target_proteins[prot_seq], sep = "\n"), x = "Residue", y = "No. of Hits") +
        			ylim(0, max = max(as.numeric(eluate_tibble[,"Repeat_No"])))
        	plot(plot)
        	plot_list <- list.append(plot_list, plot)
        	assign(toString(make.names(paste(toString(target_proteins[prot_seq]), "_graph", sep = ""))), plot)
    	}
    	
    	new_table_list <- list()
    	
    	for (graph_tables in 1:length(table_list)) {
        	new_table_list <- list.append(new_table_list, as.name(table_list[graph_tables]))
    	}
    	
    	R_target_proteins <- make.names(getName(protein_fasta))
    	
    	names(plot_list) <- R_target_proteins
    	
    	# Return the .fasta file with coloured text and appropriate legend ----
    	
    	color_pallette <- brewer.pal(max(as.numeric(eluate_tibble[,"Repeat_No"])), "Reds")
    	
    	color_pallette_key <- c()
    	
    	# Make the legend:
    	for (colors in 1:length(color_pallette)) {
    	    color_pallette_key <- c(color_pallette_key, paste("<span style=\"color:", toString(color_pallette[colors]), "\">", toString(colors), " Hit(s)", "</span>", sep = ""))
    	}
    	
    	color_pallette_key_list <- c()
    	
    	for (a in 1:length(color_pallette_key)) {
    	    color_key_part <- paste(color_pallette_key[a], "<br>", sep = "")
    	    color_pallette_key_list <- c(color_pallette_key_list, color_key_part)
    	}
    	
    	color_pallette_key_list <- c("Colour legend:<br>", color_pallette_key_list)
    	
    	for (graph_tables in 1:length(table_list)) {
    	    color_tibble <- add_column(get(table_list[[graph_tables]]), residue_color = "")
    	    for (row in 1:nrow(color_tibble)) {
    	        # If detected in 0 repeats:
    	        if (color_tibble[row, "freq"] == 0) {
    	            color_tibble[row, "residue_color"] <- toString(color_tibble[row, "residue"])
    	        } else {
    	            # Assign the number to appropriate colour:
    	            color <- color_pallette[as.numeric(color_tibble[row, "freq"])]
    	            color_tibble[row, "residue_color"] <- paste("<span style=\"color:", toString(color), "\">", toString(color_tibble[row, "residue"]), "</span>", sep = "")
    	        }
    	    }
    	    assign(toString(make.names(paste(toString(target_proteins[graph_tables]), "_graph_tibble", sep = ""))), color_tibble)
    	}
    	
    	color_string_list <- c()
    	
    	for (graph_tables in 1:length(table_list)) {
    	    color_string <- as.vector(get(table_list[[graph_tables]])$residue_color)
    	    color_string <- paste(color_string, collapse = "")
    	    assign(toString(make.names(paste(toString(target_proteins[graph_tables]), "_color_string", sep = ""))), color_string)
    	    color_string_list <- c(color_string_list, toString(make.names(paste(toString(target_proteins[graph_tables]), "_color_string", sep = ""))))
    	}
    	
    	for (strings in 1:length(color_string_list)) {
    	    color_string_list[strings] <- paste("> ", target_proteins[strings], "<br>", "<p style='word-wrap: break-word'>", get(color_string_list[strings]), "</p>", sep = "")
    	}
    	
    	names(color_string_list) <- R_target_proteins
    	
    	# create a pseudo .fasta file which includes colors for residues:
    	
    	color_fasta <- c()
    	
    	for (color_strings in 1:length(color_string_list)) {
    	    color_fasta_part <- paste(color_string_list[color_strings], "<br>", "<br>", sep = "")
    	    color_fasta <- c(color_fasta, color_fasta_part)
    	}
    	
    	# Render the plots in Shiny ----
	
    	output$repeats_plot <- {(
    	    renderPlot(get(paste(make.names(input$fasta_select), "_graph", sep = "")))
    	)}
    	
    	output$fasta_color_legend <- {(
    	        renderUI(HTML(paste(color_pallette_key_list, sep = "")))
    	)}
	
    	output$fasta_color <- {(
    	    renderUI(
    	        HTML(color_string_list[make.names(input$fasta_select)])
    	    )
    	)}
    })
}

# Run the application 
shinyApp(ui = ui, server = server)


# END OF APPLICATION