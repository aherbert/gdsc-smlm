package gdsc.smlm.data.config;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;

public class TestConfigTest
{
	@Test
	public void canChangeFieldNameAndDeserialise() throws IOException
	{
		int value = 35;
		// These two messages have the same type for field 1 but a different name
		TestConfig.Message1.Builder b1 = TestConfig.Message1.newBuilder();
		TestConfig.Message2.Builder b2 = TestConfig.Message2.newBuilder();
		b1.setField1(value);
		b2.setEntry1(value);

		TestConfig.Message1 m1 = b1.build();
		TestConfig.Message2 m2 = b2.build();

		// Write as message 1
		ByteArrayOutputStream output = new ByteArrayOutputStream();
		m1.writeDelimitedTo(output);
		output.close();
		byte[] buf = output.toByteArray();

		// Read as message 2
		ByteArrayInputStream input = new ByteArrayInputStream(buf);
		TestConfig.Message2 o2 = TestConfig.Message2.parseDelimitedFrom(input);

		// They should match as the field value is the same
		Assert.assertEquals(m2, o2);
	}

	@Test
	public void cannotChangeFieldTypeAndDeserialise() throws IOException
	{
		int value = 35;
		// These two messages have a different type (and byte size) for field 1
		TestConfig.Message1.Builder b1 = TestConfig.Message1.newBuilder();
		TestConfig.Message3.Builder b3 = TestConfig.Message3.newBuilder();
		b1.setField1(value);
		//b1.setField2(true);
		b3.setField1(value);
		//b3.setField2(1);

		TestConfig.Message1 m1 = b1.build();
		TestConfig.Message3 m3 = b3.build();

		// Write as message 1
		ByteArrayOutputStream output = new ByteArrayOutputStream();
		m1.writeDelimitedTo(output);
		output.close();
		byte[] buf = output.toByteArray();

		// Read as message 3
		ByteArrayInputStream input = new ByteArrayInputStream(buf);
		TestConfig.Message3 o3 = TestConfig.Message3.parseDelimitedFrom(input);

		// They should not match as the field size has changed
		Assert.assertNotEquals(m3, o3);
	}
}
