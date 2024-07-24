#!/usr/bin/awk -f

BEGIN {
    FS=" "; OFS="\t"; # Set input field separator to space and output field separator to tab
    print "Accession", "Location", "Completeness", "Topology", "Organism", "tissue_type", "Specimen_Voucher", "Isolate", "Gcode", "Country", "Lat_lon", "Collection_date", "Collected_by", "Identified_by", "Tech"; # Print header
}

{
    # Extract GenBank accession number
    #match($1, /gb\|([^|]+)\|/, arr);
	#match($1, /(gb|ref)\|([^|]+)\|/, arr);
	
	    #accession = arr[1];
	# Extract GenBank accession number between '>' and '|'
    match($1, />\w*\|\w*\.?\w*\|/, arr);
    gsub(/>|\|/, "", arr[0]);
    accession = arr[0];

    
    # Initialize variables for possible fields
    location = completeness = topology = organism = tissue_type = specimen_voucher = isolate = gcode = country = lat_lon = collection_date = collected_by = identified_by = tech = "";
    
    # Loop through fields starting from the second to extract key-value pairs
    for (i = 2; i <= NF; i++) {
        if ($i ~ /^\[/) {
            gsub(/\[|\]/, "", $i); # Remove brackets
            split($i, kv, "="); # Split by '=' to get key-value pair
            key = kv[1]; value = kv[2];
            
            # Assign value to corresponding variable based on key
            if (key == "location") location = value;
			else if (key == "completeness") completeness = value;
			else if (key == "topology") topology = value;
            else if (key == "organism") organism = value;
			else if (key == "tissue_type") tissue_type = value;
			else if (key == "specimen_voucher") specimen_voucher = value;
            else if (key == "isolate") isolate = value;
            else if (key == "gcode") gcode = value;
            else if (key == "country") country = value;
            else if (key == "lat_lon") lat_lon = value;
            else if (key == "collection_date") collection_date = value;
            else if (key == "collected_by") collected_by = value;
			else if (key == "identified_by") identified_by = value;
            else if (key == "tech") tech = value;
        }
    }
    
    # Print the extracted data
    print accession, location, completeness, topology, organism, tissue_type, specimen_voucher, isolate, gcode, country, lat_lon, collection_date, collected_by, identified_by, tech;
}
