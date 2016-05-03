/**
 * Created by robertopreste on 29/04/2016.
 */

function createBarplot(sitename) {
    // Create a barplot based on the properties found in the scraped website

    var margin = {top: 30, right: 40, bottom: 160, left: 50},
        width = 800 - margin.left - margin.right,
        height = 600 - margin.top - margin.bottom;

    var x = d3.scale.ordinal()
        .rangeRoundBands([0, width], .1);

    var y = d3.scale.linear()
        .range([height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var svg = d3.select("#barchart")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .attr("preserveAspectRatio", "xMinYMin meet")
        .append("g")
        .attr("transform", "translate(" + margin.left + ", " + margin.top + ")");

    d3.csv(sitename + "/tagsFound.csv", type, function (error, data) {
    //d3.csv("tagsFound.csv", type, function (error, data) {

        x.domain(data.map(function (d) { return d["Type"]; }));
        y.domain([0, d3.max(data, function (d) { return d["Value"]; })]);

        var bar = svg.selectAll("g")
            .data(data)
            .enter().append("g")
            .attr("transform", function (d, i) {
                return "translate(" + x(d["Type"]) + ", 0)";
            });

        bar.append("rect")
            .attr("y", function (d) {
                return y(d["Value"]);
            })
            .attr("height", function (d) {
                return height - y(d["Value"]);
            })
            .attr("width", x.rangeBand())
            .attr("fill", "steelblue")
            .on("mouseover", function (d) {
                d3.select(this)
                    .attr("fill", "brown")
                    .append("title")
                    .text(function (d) { return d["Type"]; })
            })
            .on("mouseout", function (d) {
                d3.select(this)
                    .attr("fill", "steelblue");
            })
            .on("click", function (d) {
                d3.select("#bubblechart")
                    .selectAll(".bubble")
                    .remove();
                bubbleMaker(d["Type"]);
                d3.select("#bubbleinfo")
                    .selectAll("table")
                    .remove();
                infoMaker(d["Type"]);
            });

        bar.append("text")
            .attr("x", x.rangeBand() / 2)
            .attr("y", function (d) {
                return y(d["Value"]) + 10;
            })
            .attr("dy", ".75em")
            .text(function (d) {
                return d["Value"];
            })
            .attr("font-family", "sans-serif")
            .attr("font-size", "16px")
            .attr("fill", "white")
            .attr("text-anchor", "middle");

        svg.append("g")
            .attr("transform", "translate(0, " + height + ")")
            .call(xAxis)
            .selectAll("text")
            .style("text-anchor", "end")
            .attr("font-family", "sans-serif")
            .attr("font-size", "14px")
            .attr("dx", "-.8em")
            .attr("dy", ".15em")
            .attr("transform", "rotate(-65)");

    });
}

function type(d) {
    d["Value"] = +d["Value"];
    return d;
}

function bubbleMaker(csvfile) {
    // Create a bubblechart based on the specific property type selected in the barplot
    
    var diameter = 400,
        color = d3.scale.category20c();

    var bubble = d3.layout.pack()
        .sort(null)
        .size([diameter, diameter])
        .padding(1.5);

    var svg = d3.select("#bubblechart")
        .append("svg")
        .attr("width", diameter)
        .attr("height", diameter)
        .attr("class", "bubble")
        .attr("id", csvfile)
        .style("padding", "30px")
        .style("margin", "10px")
        .style("background-color", "#fff");

    d3.select("#" + csvfile)
        .append("text")
        .text(csvfile)
        .attr("x", diameter / 3)
        .attr("y", -10)
        .style("font-weight", "bold")
        .style("font-size", "16px")
        .style("text-align", "center");

    d3.csv(csvfile + ".csv", function (error, data) {
        data = data.map(function (d) {
            d.value = +d["Value"];
            return d;
        });

        var nodes = bubble.nodes({children:data})
            .filter(function (d) {
                return !d.children;
            });

        var bubbles = svg.append("g")
            .attr("transform", "translate(0,0)")
            .selectAll(".bubble")
            .data(nodes)
            .enter();

        bubbles.append("circle")
            .attr("r", function (d) {
                return d.r;
            })
            .attr("cx", function (d) {
                return d.x;
            })
            .attr("cy", function (d) {
                return d.y;
            })
            .style("fill", function (d) {
                return color(d.value);
            });

        bubbles.insert("text")
            .attr("x", function (d) {
                return d.x;
            })
            .attr("y", function (d) {
                return d.y + 5;
            })
            .attr("text-anchor", "middle")
            .text(function (d) {
                return d["Property"] + " (" + d["Value"] + ")";
            })
            .style({"fill": "black",
                    "font-family": "Helvetica Neue, Helvetica, Arial, sans-serif",
                    "font-size": "12px",
                    "font-weight": "bold"});
    })
}

function infoMaker(csvfile) {
    // Create a table that shows all the details about the properties shown in the bubblechart

    var table = d3.select("#bubbleinfo")
        .append("table");
    var thead = table.append("thead");
    var tbody = table.append("tbody");
    var columns = ["Property", "Value", "Data"];

    thead.append("tr")
        .selectAll("th")
        .data(columns)
        .enter()
        .append("th")
        .text(function (columns) {return columns;})
        .style({"font-family": "Helvetica Neue, Helvetica, Arial, sans-serif",
                "font-size": "18px"
        });

    d3.csv(csvfile + ".csv", function (error, data) {
        data = data.map(function (d) {
            d.value = +d["Value"];
            return d;
        });

        var rows = tbody.selectAll("tr")
            .data(data)
            .enter()
            .append("tr");


        rows.selectAll("td")
            .data(function (row) {
                return columns.map(function (column) {
                    return {column: column, value: row[column]};
                });
            })
            .enter()
            .append("td")
            .style("font-family", "Helvetica Neue, Helvetica, Arial, sans-serif")
            .html(function (d) {
                return d.value;
            })
    });

    return table;

}