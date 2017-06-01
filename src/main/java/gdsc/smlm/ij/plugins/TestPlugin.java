package gdsc.smlm.ij.plugins;

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
public class TestPlugin implements PlugIn
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog("Test");
		gd.addChoice("Select", new String[] { "One", "Two" }, "");
		final Choice c2 = gd.addAndGetChoice("Select", new String[] { "Three", "Four" }, "");
		gd.addAndGetButton("Options", new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				ExtendedGenericDialog gd2 = new ExtendedGenericDialog("Test2", null); // This makes it model
				gd2.addMessage(c2.getSelectedItem());
				gd2.showDialog();
			}
		});
		gd.addStringField("Another", "field");
		gd.addStringField("Testing", "Some text", 15, new OptionListener<TextField>()
		{
			@Override
			public void collectOptions(TextField field)
			{
				IJ.log(field.getText());
			}
		});
		gd.addFilenameField("File", "", 30);
		gd.addChoice("Select", new String[] { "Five", "Six" }, "", new OptionListener<Choice>()
		{
			@Override
			public void collectOptions(Choice field)
			{
				IJ.log(field.getSelectedItem());
			}
		});
		gd.showDialog();
	}
}
