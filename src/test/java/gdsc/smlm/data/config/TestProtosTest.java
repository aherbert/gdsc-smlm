package gdsc.smlm.data.config;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;

public class TestProtosTest
{
	@Test
	public void canChangeFieldNameAndDeserialise() throws IOException
	{
		int value = 35;
		// These two messages have the same type for field 1 but a different name
		TestProtos.Message1.Builder b1 = TestProtos.Message1.newBuilder();
		TestProtos.Message2.Builder b2 = TestProtos.Message2.newBuilder();
		b1.setField1(value);
		b2.setEntry1(value);

		TestProtos.Message1 m1 = b1.build();
		TestProtos.Message2 m2 = b2.build();

		// Write as message 1
		ByteArrayOutputStream output = new ByteArrayOutputStream();
		m1.writeDelimitedTo(output);
		output.close();
		byte[] buf = output.toByteArray();

		// Read as message 2
		ByteArrayInputStream input = new ByteArrayInputStream(buf);
		TestProtos.Message2 o2 = TestProtos.Message2.parseDelimitedFrom(input);

		// They should match as the field value is the same
		Assert.assertEquals(m2, o2);
	}

	@Test
	public void cannotChangeFieldTypeAndDeserialise() throws IOException
	{
		int value = 35;
		// These two messages have a different type (and byte size) for field 1
		TestProtos.Message1.Builder b1 = TestProtos.Message1.newBuilder();
		TestProtos.Message3.Builder b3 = TestProtos.Message3.newBuilder();
		b1.setField1(value);
		//b1.setField2(true);
		b3.setField1(value);
		//b3.setField2(1);

		TestProtos.Message1 m1 = b1.build();
		TestProtos.Message3 m3 = b3.build();

		// Write as message 1
		ByteArrayOutputStream output = new ByteArrayOutputStream();
		m1.writeDelimitedTo(output);
		output.close();
		byte[] buf = output.toByteArray();

		// Read as message 3
		ByteArrayInputStream input = new ByteArrayInputStream(buf);
		TestProtos.Message3 o3 = TestProtos.Message3.parseDelimitedFrom(input);

		// They should not match as the field size has changed
		Assert.assertNotEquals(m3, o3);
	}
}
