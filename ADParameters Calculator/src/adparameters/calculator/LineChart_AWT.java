/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package adparameters.calculator;

import java.util.ArrayList;
import java.util.Iterator;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author tasos
 */
public class LineChart_AWT {

    private JFreeChart lineChart;

    public LineChart_AWT(String chartTitle, String xTitle, String yTitle, ArrayList xData, ArrayList yData) {
        this.lineChart = ChartFactory.createXYLineChart(
                chartTitle,
                xTitle,
                yTitle,
                createDataset(xData, yData),
                PlotOrientation.VERTICAL,
                true, true, false);
    }

    private XYSeriesCollection createDataset(ArrayList xData, ArrayList yData) {
        XYSeries series = new XYSeries("2016");
        Iterator ity = yData.iterator();
        Iterator itx = xData.iterator();
        Double x, y;
        while (ity.hasNext()) {
            x = (Double) itx.next();
            y = (Double) ity.next();
            series.add(x, y);
        }
        
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);

        return dataset;
    }

    public JFreeChart getLineChart() {
        return lineChart;
    }

}
