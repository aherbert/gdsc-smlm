/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.utils;

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
