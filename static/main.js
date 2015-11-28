$(function() {
  console.log('jquery is working!Jacot');
  createGraph();
});

function createGraph() {
  var margin = 120;
  var w = 700 - 2 * margin, h = 500 - 2 * margin;
  var svg = d3.select("#chart")
                .append("svg")
                    .attr("width", w + 2 * margin)
                    .attr("height", h + 2 * margin)
                .append("svg:g")
                        .attr("transform", "translate(" + 3*margin + ", " + margin + ")");

  var x_axis_scale = d3.scale.ordinal().rangeBands([0, w])
                            
  var xscale = d3.scale.linear().range([0, w]);
  var yscale = d3.scale.linear().range([h, 0]);

  var callback = function(data){
    var matrix = []
    var pairNames = []
    var geneNames = []
    var colSortOrderDesc;
    
    colSortOrderDesc = false;
    for (var pairName in data) {
        if (data.hasOwnProperty(pairName)) {
            pairNames.push(pairName);
            var genes= data[pairName]
            var n = genes.length
            curr_pair = []
            genes.forEach(function(gene){
                curr_pair.push(gene);
                if (geneNames.length < n){
                    geneNames.push(gene[2]);
                }
            });          
            matrix.push(curr_pair)
        }
    }   
    x_axis_scale.domain(d3.range(matrix[0].length)) // set x_axis_scale's domain to be number of genes in a pair
    var max_pval = d3.max(matrix, function(pair) {
                   return d3.max(pair.map(function(gene){return parseFloat(gene[3])})) 
                   }); 
    var min_pval = d3.min(matrix, function(pair) {
                   return d3.min(pair.map(function(gene){return parseFloat(gene[3])}));
                   }); 

    var color_scale=d3.scale.pow()
                        .exponent(0.25)
                        .domain([min_pval,max_pval])
                        .range(['red','blue']);

    svg.append("rect")
      .attr("class", "background")
      .attr("width", w)
      .attr("height", h);

    var pair = svg.selectAll(".pair")
      .data(matrix)
    .enter().append("g")
      .attr("class", "pair")
      .attr("transform", function(d, i) { return "translate(0," + x_axis_scale(i) + ")"; })
      .each(pair);

    pair.append("line")
        .attr("x2",w)
        .style("stroke","#fff");

    pair.append("text")
        .attr("x",-6)
        //.text(function(d,i) { return "haha" });
        .attr("y",x_axis_scale.rangeBand()/2)
        .attr("dy", ".32em")
        .attr("text-anchor", "end")
        .text(function(d,i) { return pairNames[i]; });

    svg.selectAll(".column").data([]).exit().remove();
    var column =svg.selectAll(".column")
              .data(geneNames)
              .enter().append("g")
              .attr("class","column")
              .attr("transform", function(d, i) { return "translate(" + i*15 + ")rotate(-90)"; });
               
    column.append("line")
  
    function getColumnI(i){
      var index=0;
      var log2r=[];
      var cd;
      for(index=0; index<matrix.length; index++){
           genes=matrix[index];
           if(genes.length>i){
                log2r.push(parseFloat(genes[i][3]));
           }else{
                log2r.push(0);
           }
      }
      return log2r;
    }

    function sortbylabel(i,colSortOrderDesc){
       var t = svg.transition().duration(3000);
       var log2r=getColumnI(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(matrix.length).sort(function(a,b){
        if(colSortOrderDesc == true){
            return log2r[a] - log2r[b];
        }else{
            return log2r[b] - log2r[a];
        }
       });
       t.selectAll(".pair")
        .attr("transform", function(d,i) { return "translate(0,"+sorted.indexOf(i) * x_axis_scale.rangeBand()+")";});
  }  


    column.append("line")
        .attr("x1",-w)
        .style("stroke","#fff");


    column.append("text")
      .attr("x", 6)
      .attr("y", x_axis_scale.rangeBand()/2)
      .attr("width",x_axis_scale.rangeBand())
      .attr("height",x_axis_scale.rangeBand())
      .attr("dy", ".32em")
      .text(function(d, i) { return geneNames[i]; }) 
      .on("click", function(d,i) { colSortOrderDesc=!colSortOrderDesc;sortbylabel(i,colSortOrderDesc);}) ;


 
  function pair(pair) {
    var cell = d3.select(this).selectAll(".gene")
        .data(pair)
      .enter().append("rect")
        .attr("class", "gene")
        .attr("x", function(d,i) { return x_axis_scale(i); })
        .attr("width", x_axis_scale.rangeBand())
        .attr("height", x_axis_scale.rangeBand())
        .attr("title",function(d){return d[3]})
        .style("fill", function(d) { return color_scale(parseFloat(d[3])) });
    }

    console.log('you called callback! you know how to get data!');
  };

  d3.json("/data", callback);

  // Code goes here
}
