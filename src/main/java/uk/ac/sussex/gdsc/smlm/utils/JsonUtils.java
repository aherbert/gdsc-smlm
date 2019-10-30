/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.utils;

/**
 * Provides utilities for working with JSON.
 */
public final class JsonUtils {
  /** Allowed characters. */
  private static final char[] ALLOWED = {'-', '_'};

  /** No public constructor. */
  private JsonUtils() {}

  /**
   * Simplify the JSON string. Removes redundant double quotes around fields, e.g. if they have only
   * letters or digits.
   *
   * <p>If the input is null then an empty string will be returned.
   *
   * @param json the json
   * @return the simplified string
   */
  public static String simplify(String json) {
    if (json == null || json.length() == 0) {
      return "";
    }
    final char[] chars = json.toCharArray();
    final char[] newChars = new char[chars.length];
    int length = 0;
    // Scan for first double quote
    int index = 0;
    final int size = chars.length;
    while (index < size) {
      if (isDoubleQuote(chars, index)) {
        final int start = index;
        // Scan for end double quote
        int end = -1;
        for (int j = index + 1; j < size; j++) {
          if (isDoubleQuote(chars, j)) {
            end = j;
            break;
          }
        }
        if (end > start) {
          // We have a terminating double quote.
          // Check if all the characters between start and end are valid
          if (canSimplify(chars, start + 1, end)) {
            for (index = start + 1; index < end; index++) {
              newChars[length++] = chars[index];
            }
            index = end + 1;
            continue;
          }
          // Cannot remove the quotes so advance to copy all the chars
          end++;
        } else {
          // No terminating double quote so this is the end of the string
          end = size;
        }

        // We cannot simplify so copy all the characters
        for (index = start; index < end; index++) {
          newChars[length++] = chars[index];
        }
      } else {
        newChars[length++] = chars[index++];
      }
    }
    return new String(newChars, 0, length);
  }

  /**
   * Checks if the character at the index is a double-quote that is not escaped (prefixed with
   * {@code '\\'}).
   *
   * @param chars the chars
   * @param index the index
   * @return true, if is double quote
   */
  private static boolean isDoubleQuote(char[] chars, int index) {
    return chars[index] == '"' && (index == 0 || chars[index - 1] != '\\');
  }

  /**
   * Check if the range of characters can be simplified, i.e. the outer double-quotes around the
   * characters can be removed.
   *
   * @param chars the characters
   * @param start the start (inclusive)
   * @param end the end (exclusive)
   * @return true if double-quotes are redundant around the characters
   */
  private static boolean canSimplify(char[] chars, int start, int end) {
    for (int j = start; j < end; j++) {
      final char ch = chars[j];
      if (Character.isWhitespace(ch)) {
        return false;
      }
      if (Character.isLetterOrDigit(ch) || isAllowed(ch)) {
        continue;
      }
      return false;
    }
    return true;
  }

  /**
   * Checks if the character is allowed.
   *
   * @param ch the character
   * @return true if allowed
   */
  private static boolean isAllowed(char ch) {
    for (int i = 0; i < ALLOWED.length; i++) {
      if (ch == ALLOWED[i]) {
        return true;
      }
    }
    return false;
  }

  /**
   * Formats a JSON string for pretty printing.
   *
   * <p>Copied from <a href="https://stackoverflow.com/a/49564514">StackOverflow pretty-print JSON
   * in Java</a>.
   *
   * <Blockquote><p>"The basic idea is to trigger the formatting based on special characters in
   * JSON. For example, if a '{' or '[' is observed, the code will create a new line and increase
   * the indent level.</p>
   *
   * <p>Disclaimer: I only tested this for some simple JSON cases (basic key-value pair, list,
   * nested JSON) so it may need some work for more general JSON text, like string value with quotes
   * inside, or special characters (\n, \t etc.)."</p></blockquote>
   *
   * @param json the JSON
   * @return pretty-print formatted JSON
   */
  public static String prettyPrintJson(String json) {
    final StringBuilder sb = new StringBuilder();
    int indentLevel = 0;
    boolean inQuote = false;
    for (char ch : json.toCharArray()) {
      switch (ch) {
        case '"':
          // switch the quoting status
          inQuote = !inQuote;
          sb.append(ch);
          break;
        case ' ':
          // For space: ignore the space if it is not being quoted.
          if (inQuote) {
            sb.append(ch);
          }
          break;
        case '{':
        case '[':
          // Starting a new block: increase the indent level
          sb.append(ch);
          indentLevel++;
          appendIndentedNewLine(indentLevel, sb);
          break;
        case '}':
        case ']':
          // Ending a new block; decrease the indent level
          indentLevel--;
          appendIndentedNewLine(indentLevel, sb);
          sb.append(ch);
          break;
        case ',':
          // Ending a json item; create a new line after
          sb.append(ch);
          if (!inQuote) {
            appendIndentedNewLine(indentLevel, sb);
          }
          break;
        default:
          sb.append(ch);
      }
    }
    return sb.toString();
  }

  /**
   * Print a new line with indentation at the beginning of the new line.
   *
   * @param indentLevel the indent level
   * @param sb the string builder
   */
  private static void appendIndentedNewLine(int indentLevel, StringBuilder sb) {
    sb.append('\n');
    for (int i = 0; i < indentLevel; i++) {
      // Assuming indentation using 2 spaces
      sb.append("  ");
    }
  }
}
