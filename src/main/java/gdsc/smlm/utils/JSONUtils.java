package gdsc.smlm.utils;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Provides utilities for working with JSON
 */
public class JSONUtils
{
	/**
	 * Simplify the JSON string. Removes redundant double quotes around fields, e.g. if they have only letters or
	 * digits.
	 * <p>
	 * If the input is null then an empty string will be returned.
	 *
	 * @param json
	 *            the json
	 * @return the simplified string
	 */
	public static String simplify(String json)
	{
		if (json == null || json.length() == 0)
			return "";
		char[] chars = json.toCharArray();
		char[] newChars = new char[chars.length];
		int length = 0;
		// Scan for first double quote
		int i = 0, size = chars.length;
		while (i < size)
		{
			if (isDoubleQuote(chars, i))
			{
				int start = i;
				// Scan for end double quote
				int end = -1;
				for (int j = i + 1; j < size; j++)
				{
					if (isDoubleQuote(chars, j))
					{
						end = j;
						break;
					}
				}
				if (end > start)
				{
					// We have a terminating double quote.
					// Check if all the characters between start and end are valid
					if (canSimplify(chars, start + 1, end))
					{
						for (i = start + 1; i < end; i++)
							newChars[length++] = chars[i];
						i = end + 1;
						continue;
					}
					else
					{
						// Cannot remove the quotes so advance to copy all the chars
						end++;
					}
				}
				else
				{
					// No terminating double quote so this is the end of the string
					end = size;
				}

				// We cannot simplify so copy all the characters
				for (i = start; i < end; i++)
					newChars[length++] = chars[i];
			}
			else
			{
				newChars[length++] = chars[i++];
			}
		}
		return new String(newChars, 0, length);
	}

	private static boolean isDoubleQuote(char[] chars, int i)
	{
		// Check the quote is not escaped
		return chars[i] == '"' && (i == 0 || chars[i - 1] != '\\');
	}

	private static boolean canSimplify(char[] chars, int start, int end)
	{
		for (int j = start; j < end; j++)
		{
			if (!Character.isLetterOrDigit(chars[j]))
				return false;
		}
		return true;
	}
}
