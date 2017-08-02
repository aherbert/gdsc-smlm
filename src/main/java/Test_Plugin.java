
import java.awt.Choice;
import java.awt.TextField;
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
		gd.addStringField("Testing", textFields[1], 15, new OptionListener<TextField>()
		{
			public boolean collectOptions(TextField field)
			{
				IJ.log(field.getText());
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
		gd.addChoice("Select3", new String[] { "Five", "Six" }, optionFields[2], new OptionListener<Choice>()
		{
			public boolean collectOptions(Choice field)
			{
				IJ.log(field.getSelectedItem());
				return true;
			}

			public boolean collectOptions()
			{
				IJ.log(optionFields[2]);
				return true;
			}
		});
		gd.showDialog();
		optionFields[0] = gd.getNextChoice();
		optionFields[1] = gd.getNextChoice();
		textFields[0] = gd.getNextString();
		textFields[1] = gd.getNextString();
		optionFields[2] = gd.getNextChoice();
		gd.collectOptions();
	}
}
