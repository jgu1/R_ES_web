$(function() {
  createGraph();
});

function createGraph() {
  
  //var x_axis_scale = d3.scale.ordinal().rangeBands([0, w])
 
  //var xscale = d3.scale.linear().range([0, w]);
  //var yscale = d3.scale.linear().range([h, 0]);


  var sub_clusterscallback = function(data){
    var gene_p_qs = data.gene_p_qs
    var pairname_idx = data.pairname_idx
   
    for (var i = 0; i < pairname_idx.length; i++){
        var curr_rowcomb_cols = pairname_idx[i];
        var row_comb = curr_rowcomb_cols[0];
        var cols = curr_rowcomb_cols[1];
        
        var cluster = {};
        // loop through all rows
        for (var i_pairname = 0; i_pairname < row_comb.length; i_pairname++){
            var curr_pairname = row_comb[i_pairname];
            var curr_all_genes = gene_p_qs[curr_pairname];
            curr_pair_genes_arr = new Array();
            // loop though all columns
            for (var i_pass_idx = 0; i_pass_idx < cols.length; i_pass_idx++){
                var curr_pass_idx = cols[i_pass_idx];
                var curr_pass_gene = curr_all_genes[curr_pass_idx];
                curr_pair_genes_arr.push(curr_pass_gene);
            }
            cluster[curr_pairname] = curr_pair_genes_arr;
        }   
        
        var a = 1;
    }
    
    function 


 
    pairname_idx.forEach(function (row_comb_cols,i){
        var cluster = new Array();
        var row_comb = row_comb_cols[0];
        var cols = row_comb_cols[1];
        

    });
    

    var a = 1
  }

  var detailcallback = function(data){
    var matrix = [];
    var pairNames = [];
    var size = 15; 
    var gene = data.gene;
    var all_SNPs_list = data.all_SNPs_list
    var pair_SNP_dict = data.pair_SNP_dict;
    var curr_pair = [];

    for (var pairName in pair_SNP_dict) {
        pairNames.push(pairName);
        var SNPs= pair_SNP_dict[pairName];
        var n = SNPs.length;
        curr_pair = []
        SNPs.forEach(function(SNP){
            curr_pair.push(SNP);
        });  
        matrix.push(curr_pair); 
    }
    var detail_height = matrix.length * size;
    var detail_width = matrix[0].length * size; // should not be null, because at least one gene has SNPs
    var margin = {top: 120, right: 200, bottom: 10, left: 360};
    
    d3.select("#detail").html("");
    d3.select("#detail").append("div") 
        .attr("id","gene_of_interest")
        .text("Gene of Interest: " + gene)
        .style("margin-left",margin.left +"px");

    var canvas2=d3.select("#detail").append("svg")
           .attr("width", detail_width + margin.left)
           .attr("height", detail_height + margin.top)
   
    //d3.select("#gene_of_interest")
    
    //canvas2.html("");// clear canvas2 before drawing a new one
    var svg2=canvas2
       .append("g")
       .attr("transform", "translate(" + margin.left + "," + margin.top + ")"); 

    var max_pval = 5;
    var min_pval = 0;
    var color_scale=d3.scale.linear()
                        .domain([min_pval,max_pval])
                        .range(['#ffeda0','#f03b20']);

    var height_detail = matrix.length * size;
    var width_detail= matrix[0].length * size;

    svg2.append("rect")
      .attr("class", "background")
      .attr("width", width_detail)
      .attr("height", height_detail);

    var pair = svg2.selectAll(".pair2")
      .data(matrix)
    .enter().append("g")
      .attr("class", "pair2")
      .attr("transform", function(d, i) { return "translate(0," + i*size + ")"; })
      .each(pair2);

     pair.append("line")
        .attr("x2",width_detail)
        .style("stroke","#fff");

     pair.append("text")
        .attr("x",-6)
        .attr("y",7)
        .attr("dy", ".32em")
        .attr("text-anchor", "end")
        .text(function(d,i) { return pairNames[i]; })
        .on("click", function(d,i){ 
            sortAlongRow(i);
        });

    svg2.selectAll(".column").data([]).exit().remove();
    var column =svg2.selectAll(".column")
              .data(all_SNPs_list)
              .enter().append("g")
              .attr("class","column")
              .attr("transform", function(d, i) { return "translate(" + i*size + ")rotate(-90)"; });
               
    column.append("line")
          .style("stroke","#fff");
 

    column.append("line")
        .attr("x1",-height_detail)
        .style("stroke","#fff");


    column.append("text")
      .attr("x", 6)
      .attr("y", 7)
      .attr("width",size)
      .attr("height",size)
      .attr("dy", ".32em")
      .text(function(d, i) { return all_SNPs_list[i]; }) 
      .on("click", function(d,i) {sortAlongColumn(i);}) ; 
    
    function pair2(pair) {
        var cell = d3.select(this).selectAll(".SNP")
            .data(pair)
          .enter().append("rect")
            .attr("class", "SNP")
            .attr("x", function(d,i) { return i*size; })
            .attr("width", size)
            .attr("height", size)
            .style("fill", 
                function(d) { 
                    if(parseFloat(d[2]) < 0){
                        return '#000'
                    }else{
                        return color_scale(-Math.log10(parseFloat(d[2]))) 
                    }
                })
            .append("title").text(function(d){
             if (parseFloat(d[2]) < 0){
                return '';   
             }else{
                return "GWAS pval: " + d[2] + "\neSNP name: " + d[1] ;
                //return "GWAS pval :" + d[2] + "\nGWAS_SNP:" + d[0] + "\neQTL_SNP :" + d[1] ;
             }
            });  
    }
 
    function getRowINegLog10(i){
      var index=0;
      var negLog10=[];
      var rowIRaw = matrix[i];
      for (var i=0; i<rowIRaw.length; i++){
        SNP = rowIRaw[i];
        if (parseFloat(SNP[2]) < 0){
            negLog10.push(-1);    
        }else{
            negLog10.push(-Math.log10(parseFloat(SNP[2])));  
        }
      }
      return negLog10;
    }


    function sortAlongRow(i){
       var t = svg2.transition().duration(3000);
       var rowIData=getRowINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(rowIData.length).sort(function(a,b){   //index array pointing to real data array
            return rowIData[b] - rowIData[a];
       });
       // sort column names
       t.selectAll(".column")
        .attr("transform", function(d,i) { return "translate(" + sorted.indexOf(i) * size+")rotate(-90)";});
       // sort rects in each row 
       var SNPs = t.selectAll(".pair2")
                    .selectAll(".SNP")
                    .attr("x", function(d,i) { return sorted.indexOf(i)*size; });
    }  

    function getColumnINegLog10(i){
      var index=0;
      var negLog10=[];
      var cd;
      for(index=0; index<matrix.length; index++){
           SNPs=matrix[index];
           if(SNPs.length>i){
                if (parseFloat(SNPs[i][2]) < 0){
                    negLog10.push(-1);
                }else{
                     negLog10.push(-Math.log10(parseFloat(SNPs[i][2])));
                }
           }else{
                negLog10.push(0);
           }
      }
      return negLog10;
    }

    function sortAlongColumn(i){
       var t = svg2.transition().duration(3000);
       var ColumnIData=getColumnINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(matrix.length).sort(function(a,b){   //index array pointing to real data arra
            return ColumnIData[b] - ColumnIData[a];
       });
       t.selectAll(".pair2")
        .attr("transform", function(d,i) { return "translate(0,"+sorted.indexOf(i) * size+")";});
  }  
  }

  var draw_pair = function(draw_pair_json_obj){
    var matrix = [];
    var pairNames = [];
    var geneNames = [];
    var size = 15;
    var geneNames = draw_pair_json_obj.filtered_gene_names; 
    var gene_p_qs = draw_pair_json_obj.gene_p_qs; 
    var pairNames = draw_pair_json_obj.sorted_pair_names; 
    for (var i = 0; i<pairNames.length;i++) {
        pairName = pairNames[i];
        var genes= gene_p_qs[pairName]
        var n = genes.length
        curr_pair = []
        genes.forEach(function(gene){
            curr_pair.push(gene);
        });          
        matrix.push(curr_pair)
    
    }   
    var max_pval = 5;
    var min_pval = 0;
    var color_scale=d3.scale.linear()
                        .domain([min_pval,max_pval])
                        .range(['#ffeda0','#f03b20']);


    var height_pair = matrix.length * size;
    var width_pair= matrix[0].length * size;
    var margin = 120;
    
    var chart = d3.select("#chart");
    //var chartHeight = chart.style("height");
    //var chartHeight = parseInt(chartHeight.substring(0,chartHeight.length-2));
    if (height_pair +  margin > 500){
        chart.attr("height",500);
    }else{
        chart.attr("height", height_pair + margin);
    }
  
    //var w = width_pair - 2 * margin, h = height_pair - 2 * margin;
    var w = width_pair, h = height_pair;
    var chart = d3.select("#chart");
    chart.selectAll("svg").remove(); 
    //var svg = d3.select("#chart")
    var svg = chart.append("svg")
                    .attr("width", w + 3 * margin)
                    .attr("height",height_pair + margin)
                    .style("overflow","scroll")
                  .append("svg:g")
                        .attr("transform", "translate(" + 3*margin + ", " + margin + ")");



    svg.append("rect")
      .attr("class", "background")
      .attr("width", width_pair)
      .attr("height", height_pair);

    chart.append("input")
         .attr("type","button")
         .attr("value","Discover sub_clusters")
         .style("margin-left",3*margin +"px")
         .on("click",function(d){
            d3.json("/sub_clusters",sub_clusterscallback);
            });


    var pair = svg.selectAll(".pair")
      .data(matrix)
    .enter().append("g")
      .attr("class", "pair")
      .attr("transform", function(d, i) { return "translate(0," + i*size + ")"; })
      .each(pair);

    pair.append("line")
        .attr("x2",w)
        .style("stroke","#fff");
 
    function getRowINegLog10(i){
      var index=0;
      var negLog10=[];
      var rowIRaw = matrix[i];
      for (var i=0; i<rowIRaw.length; i++){
        gene = rowIRaw[i];
        pval=parseFloat(gene[3]);
        if (pval < 0){
            negLog10.push(-1);
        }else{
            negLog10.push(-Math.log10(pval));  
        }  
      }
      return negLog10;
    }

    /*
    function sortEachRow(row){
        var rects = row.selectAll(".gene")
            .attr("transform",function(d,i){return "translate("+(30-i) * size+")";});
    }
    */

    function sortAlongRow(i){
       var t = svg.transition().duration(3000);
       var rowIData=getRowINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(rowIData.length).sort(function(a,b){   //index array pointing to real data array
        /*
        if(rowSortOrderDesc == true){
            return rowIData[a] - rowIData[b];
        }else{
            return rowIData[b] - rowIData[a];
        }
        */
            return rowIData[b] - rowIData[a];
       });
       // sort column names
       t.selectAll(".column")
        .attr("transform", function(d,i) { return "translate(" + sorted.indexOf(i) * size+")rotate(-90)";});
       // sort rects in each row 
       var genes = t.selectAll(".pair")
                    .selectAll(".gene")
                    .attr("x", function(d,i) { return sorted.indexOf(i)*size; });
    }  

    pair.append("text")
        .attr("x",-6)
        .attr("y",7)
        .attr("dy", ".32em")
        .attr("text-anchor", "end")
        .text(function(d,i) { return pairNames[i]; })
        .on("click", function(d,i){ 
            sortAlongRow(i);
        });
        /*
        .on("click", function(d,i) {
            d3.json("/sortAlongPair?sortAlongPairName="+pairNames[i],draw_pair);
          });
        */
    svg.selectAll(".column").data([]).exit().remove();
    var column =svg.selectAll(".column")
              .data(geneNames)
              .enter().append("g")
              .attr("class","column")
              .attr("transform", function(d, i) { return "translate(" + i*size + ")rotate(-90)"; });
               
    column.append("line")
          .style("stroke","#fff");
    
    function getColumnINegLog10(i){
      var index=0;
      var negLog10=[];
      var cd;
      for(index=0; index<matrix.length; index++){
           genes=matrix[index];
           if(genes.length>i){
                columnINumber = parseFloat(genes[i][3]);
                if (columnINumber < 0){
                    negLog10.push(0);
                }else{
                    negLog10.push(-Math.log10(columnINumber));
                }
           }else{
                negLog10.push(0);
           }
      }
      return negLog10;
    }

    function sortAlongColumn(i){
       var t = svg.transition().duration(3000);
       var ColumnIData=getColumnINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(matrix.length).sort(function(a,b){   //index array pointing to real data arra
            return ColumnIData[b] - ColumnIData[a];
       });
       t.selectAll(".pair")
        .attr("transform", function(d,i) { return "translate(0,"+sorted.indexOf(i) * size+")";});
  }  


    column.append("line")
        .attr("x1",-w)
        .style("stroke","#fff");


    column.append("text")
      .attr("x", 6)
      .attr("y", 7)
      .attr("width",size)
      .attr("height",size)
      .attr("dy", ".32em")
      .text(function(d, i) { return geneNames[i]; }) 
      .on("click", function(d,i) {sortAlongColumn(i);}) ;
 
  function pair(pair) {
    var cell = d3.select(this).selectAll(".gene")
        .data(pair)
      .enter().append("rect")
        .attr("class", "gene")
        /*
        .attr("x", function(d,i) { return x_axis_scale(i); })
        .attr("width", x_axis_scale.rangeBand())
        .attr("height", x_axis_scale.rangeBand())
        */
        .attr("x", function(d,i) { return i*size; })
        .attr("width", size)
        .attr("height", size)
        .attr("title",function(d){return d[3]})
        .style("fill", 
            function(d) { 
                if(parseFloat(d[3]) < 0){
                    return '#000'
                }else{
                    return color_scale(-Math.log10(parseFloat(d[3]))) 
                }
            })
        .on("click",function(d){
            //d3.json("/detail?GWAS="+d[0]+"&eQTL="+d[1]+"&gene="+d[2],detailcallback);
            d3.json("/detail?gene="+d[2],detailcallback);
            //plotDetail(d);
            }
           )
        .append("title")
            .text(function(d){
                return "pval:" + d[3] + "\nqval:" + d[4] + "\nGWAS:" + d[0] + "\neQTL:" + d[1] + "\ngene:" + d[2];
                
            });
     
  }

  function plotDetail(name){
    var a = 1;
  }
    console.log('you called callback! you know how to get data!');
  };

  var draw_pair_json_text = document.getElementById("draw_pair_json_obj").value;
  var draw_pair_json_obj = JSON.parse(draw_pair_json_text);
  //var filtered_gene_names_json_text = document.getElementById("filtered_gene_names_json_obj").value;
  //var filtered_gene_names_json_obj = JSON.parse(filtered_gene_names_json_text);
  //draw_pair(gene_p_qs_json_obj,filtered_gene_names_json_obj);
  draw_pair(draw_pair_json_obj);
    // Code goes here
}
