while IFS= read -r species; 
do     
species=$(echo "$species" | tr -d '\r');     
taxid=$(esearch -db taxonomy -query "$species" < /dev/null 2>/dev/null \
      | efetch -format xml 2>/dev/null \
      | xtract -pattern Taxon -element TaxId 2>/dev/null || echo "ERROR");     
echo "$species|$taxid";     
sleep 0.5; 
done < BGE_species_uniq.csv
