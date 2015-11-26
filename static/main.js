$(function() {
  console.log('jquery is working!Jacot');
  createGraph();
});

function createGraph() {
  var margin = 80;
  var w = 700 - 2 * margin, h = 500 - 2 * margin;
  var svg = d3.select("#chart")
                .append("svg")
                    .attr("width", w + 2 * margin)
                    .attr("height", h + 2 * margin)
                .append("svg:g")
                        .attr("transform", "translate(" + margin + ", " + margin + ")");

  var x_axis_scale = d3.scale.ordinal().rangeBands([0, w])
                            
  var xscale = d3.scale.linear().range([0, w]);
  var yscale = d3.scale.linear().range([h, 0]);


  svg.append("svg:g").attr("class", "x axisyo")
    .attr("transform", "translate(0, " + h + ")");
  svg.append("svg:g").attr("class", "y axis");

  var callback = function(data){
    var matrix = []
    var pairNames = []
    for (var pairName in data) {
        if (data.hasOwnProperty(pairName)) {
            pairNames.push(pairName);
            var genes= data[pairName]
            var n = genes.length
            curr_pair = []
            genes.forEach(function(gene){
                curr_pair.push(gene);
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

    pair.append("text")
        .attr("x",-80)
        //.text(function(d,i) { return "haha" });
        .attr("y",x_axis_scale.rangeBand())
        .text(function(d,i) { return pairNames[i]; });

    var gene =svg.selectAll(".gene")
        


    gene.append("text")
      .attr("x", 6)
      .attr("y", x.rangeBand() / 2)
      .attr("dy", ".32em")
      .attr("text-anchor", "start")
      .text(function(d, i) { return nodes[i].name; });


    function pair(pair) {
    var cell = d3.select(this).selectAll(".gene")
        .data(pair)
      .enter().append("rect")
        .attr("class", "gene")
        .attr("x", function(d,i) { return x_axis_scale(i); })
        .attr("width", x_axis_scale.rangeBand())
        .attr("height", x_axis_scale.rangeBand())
        .attr("title",function(d){return d[3]})
        .style("fill", function(d) { return color_scale(parseFloat(d[3])) })
           
    }



    console.log('you called callback! you know how to get data!');
  };

  d3.json("/data", callback);

  // Code goes here
}
