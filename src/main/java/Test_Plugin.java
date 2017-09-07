
import java.awt.Choice;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import ij.IJ;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ExtendedGenericDialog.OptionListener;
import ij.plugin.PlugIn;

/**
 * A simple class used to test plugin functionality
 */
public class Test_Plugin implements PlugIn
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		// The parameters that have options must be available statically for the OptionListener
		final String[] textFields = { "Some text", "More text" };
		final String[] optionFields = { "", "", "" };
		final double[] numberFields = { 2.567, 7, 4.567 };

		ExtendedGenericDialog gd = new ExtendedGenericDialog("Test");
		gd.addChoice("Select1", new String[] { "One", "Two" }, optionFields[0]);
		final Choice c2 = gd.addAndGetChoice("Select2", new String[] { "Three", "Four" }, optionFields[1]);
		gd.addAndGetButton("Options", new ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				ExtendedGenericDialog gd2 = new ExtendedGenericDialog("Test2", null); // This makes it model
				gd2.addMessage(c2.getSelectedItem());
				gd2.showDialog(true);
				gd2.getNextChoice();
			}
		});
		gd.addStringField("Another", textFields[0]);
		gd.addStringField("Testing", textFields[1], 15, new OptionListener<String>()
		{
			public boolean collectOptions(String field)
			{
				IJ.log(field);
				return true;
			}

			public boolean collectOptions()
			{
				IJ.log(textFields[1]);
				return true;
			}
		});
		gd.addFilenameField("File", "", 30);
		gd.addDirectoryField("Dir", "", 30);
		gd.addChoice("Select3", new String[] { "Five", "Six" }, optionFields[2], new OptionListener<Integer>()
		{
			public boolean collectOptions(Integer field)
			{
				IJ.log(Integer.toString(field));
				return true;
			}

			public boolean collectOptions()
			{
				IJ.log(optionFields[2]);
				return true;
			}
		});
		gd.addSlider("Slider1", 0.5, 4.5, numberFields[0], new OptionListener<Double>(){

			public boolean collectOptions(Double field)
			{
				IJ.log(field.toString());
				return true;
			}

			public boolean collectOptions()
			{
				IJ.log(Double.toString(numberFields[0]));
				return true;
			}});
		gd.addSlider("Slider2", 0, 10, numberFields[1], new OptionListener<Double>(){

			public boolean collectOptions(Double field)
			{
				IJ.log(field.toString());
				return true;
			}

			public boolean collectOptions()
			{
				IJ.log(Double.toString(numberFields[1]));
				return true;
			}});
		gd.addNumericField("Number1", numberFields[2], 2, new OptionListener<Double>(){

			public boolean collectOptions(Double field)
			{
				IJ.log(field.toString());
				return true;
			}

			public boolean collectOptions()
			{
				IJ.log(Double.toString(numberFields[2]));
				return true;
			}});
		gd.setMaxUnscrolledSize(0, 300);
		gd.showDialog();
		optionFields[0] = gd.getNextChoice();
		optionFields[1] = gd.getNextChoice();
		textFields[0] = gd.getNextString();
		textFields[1] = gd.getNextString();
		optionFields[2] = gd.getNextChoice();
		numberFields[0] = gd.getNextNumber();
		numberFields[1] = gd.getNextNumber();
		numberFields[2] = gd.getNextNumber();
		gd.collectOptions();
	}
}
